import React, { useState, useEffect, useMemo, useRef } from 'react';
import { Play, RotateCcw, Plus, Trash2, Settings, ChevronDown, ChevronUp, AlertCircle, Download } from 'lucide-react';

// --- TYPES ---

type Span = {
    id: string;
    length: number;
};

type Axle = {
    id: string;
    load: number;
    spacing: number; // Spacing to the NEXT axle
};

type AnalysisConfig = {
    E: number; // Pascal
    I: number; // m^4
    nElemsPerSpan: number;
    truckIncrement: number;
    loadCase: 'truck' | 'lane';
};

type EnvelopePoint = {
    x: number;
    max: number;
    min: number;
};

type ReactionEnvelope = {
    x: number;
    max: number;
    min: number;
};

type AnalysisResults = {
    shear: EnvelopePoint[];
    moment: EnvelopePoint[];
    deflection: EnvelopePoint[];
    xNodes: number[];
    reactions: ReactionEnvelope[];
    supportPositions: number[];
};

// --- CONSTANTS ---

const DEFAULT_SPANS: Span[] = [
    { id: 's1', length: 20 },
    { id: 's2', length: 25 },
    { id: 's3', length: 20 },
];

const DEFAULT_AXLES: Axle[] = [
    { id: 'a1', load: 50, spacing: 3.6 },
    { id: 'a2', load: 125, spacing: 1.2 },
    { id: 'a3', load: 125, spacing: 6.6 },
    { id: 'a4', load: 175, spacing: 6.6 },
    { id: 'a5', load: 150, spacing: 0 },
];

const DEFAULT_CONFIG: AnalysisConfig = {
    E: 200000000000,
    I: 0.005,
    nElemsPerSpan: 32,
    truckIncrement: 0.5,
    loadCase: 'truck',
};

// --- FEM ENGINE (Ported from VBA) ---

class BeamFEM {
    private config: AnalysisConfig;
    private spans: Span[];
    private axles: Axle[];

    constructor(spans: Span[], axles: Axle[], config: AnalysisConfig) {
        this.spans = spans;
        this.axles = axles;
        this.config = config;
    }

    // Gaussian Elimination Solver
    private solveSystem(K: number[][], F: number[], constrained: boolean[]): number[] {
        const n = F.length;
        const map: number[] = [];

        // Map free DOFs
        for (let i = 0; i < n; i++) {
            if (!constrained[i]) map.push(i);
        }

        const nFree = map.length;
        const K_red: number[][] = Array(nFree).fill(0).map(() => Array(nFree).fill(0));
        const F_red: number[] = Array(nFree).fill(0);

        // Build reduced system
        for (let i = 0; i < nFree; i++) {
            F_red[i] = F[map[i]];
            for (let j = 0; j < nFree; j++) {
                K_red[i][j] = K[map[i]][map[j]];
            }
        }

        // Forward Elimination
        for (let k = 0; k < nFree - 1; k++) {
            for (let i = k + 1; i < nFree; i++) {
                const factor = K_red[i][k] / K_red[k][k];
                for (let j = k; j < nFree; j++) {
                    K_red[i][j] -= factor * K_red[k][j];
                }
                F_red[i] -= factor * F_red[k];
            }
        }

        // Back Substitution
        const U_red: number[] = Array(nFree).fill(0);
        for (let i = nFree - 1; i >= 0; i--) {
            let sum = 0;
            for (let j = i + 1; j < nFree; j++) {
                sum += K_red[i][j] * U_red[j];
            }
            U_red[i] = (F_red[i] - sum) / K_red[i][i];
        }

        // Map back to full vector
        const U_Full: number[] = Array(n).fill(0);
        for (let i = 0; i < nFree; i++) {
            U_Full[map[i]] = U_red[i];
        }

        return U_Full;
    }

    public runAnalysis(): AnalysisResults {
        const { E, I, nElemsPerSpan, truckIncrement, loadCase } = this.config;
        const numSpans = this.spans.length;
        const nTotalElems = numSpans * nElemsPerSpan;
        const nNodes = nTotalElems + 1;
        const nDOF = nNodes * 2;

        // 1. Mesh Generation
        const elemLens: number[] = [];
        const xNodes: number[] = [0];
        let currentX = 0;

        // Stepped Shear Mesh (Start/End of each element)
        const xShear: number[] = [];

        for (const span of this.spans) {
            const le = span.length / nElemsPerSpan;
            for (let i = 0; i < nElemsPerSpan; i++) {
                xShear.push(currentX);       // Start of elem
                currentX += le;
                xShear.push(currentX);       // End of elem

                xNodes.push(currentX);
                elemLens.push(le);
            }
        }

        // 2. Global Stiffness Matrix
        const K_Global: number[][] = Array(nDOF).fill(0).map(() => Array(nDOF).fill(0));

        for (let i = 0; i < nTotalElems; i++) {
            const le = elemLens[i];
            const coeff = (E * I) / Math.pow(le, 3);
            const k_loc = [
                [12, 6 * le, -12, 6 * le],
                [6 * le, 4 * le * le, -6 * le, 2 * le * le],
                [-12, -6 * le, 12, -6 * le],
                [6 * le, 2 * le * le, -6 * le, 4 * le * le]
            ];

            const mapDOFs = [i * 2, i * 2 + 1, (i + 1) * 2, (i + 1) * 2 + 1];

            for (let r = 0; r < 4; r++) {
                for (let c = 0; c < 4; c++) {
                    K_Global[mapDOFs[r]][mapDOFs[c]] += k_loc[r][c] * coeff;
                }
            }
        }

        // 3. Boundary Conditions (Pins) and Support Positions
        const constrained: boolean[] = Array(nDOF).fill(false);
        const supportPositions: number[] = [0]; // First support at x=0
        let nodeIdx = 0;
        constrained[0] = true; // First node pinned Y
        let cumLen = 0;
        for (let i = 0; i < numSpans; i++) {
            cumLen += this.spans[i].length;
            supportPositions.push(cumLen);
            nodeIdx += nElemsPerSpan;
            constrained[nodeIdx * 2] = true; // End of span pinned Y
        }

        // 4. Initialize Envelopes
        const initEnvelope = (len: number) => Array(len).fill(0).map(() => ({ max: -Infinity, min: Infinity }));

        // Shear uses the stepped mesh (2 points per element)
        const shearEnv = initEnvelope(nTotalElems * 2).map((p, i) => ({ ...p, x: xShear[i] }));

        // Moment/Deflection use nodal mesh
        const momEnv = initEnvelope(nNodes).map((p, i) => ({ ...p, x: xNodes[i] }));
        const defEnv = initEnvelope(nNodes).map((p, i) => ({ ...p, x: xNodes[i] }));

        // Reaction envelopes (one per support)
        const reactionEnv: ReactionEnvelope[] = supportPositions.map(x => ({ x, max: -Infinity, min: Infinity }));

        // Prepare Loads
        let activeAxles = [...this.axles];
        let w_udl = 0;

        if (loadCase === 'lane') {
            w_udl = 9; // kN/m
            activeAxles = activeAxles.map(a => ({ ...a, load: a.load * 0.8 }));
        } else {
            const dla = 0.25;
            activeAxles = activeAxles.map(a => ({ ...a, load: a.load * (1 + dla) }));
        }

        // 5. UDL Analysis (Simplified Worst Case approximation for demo: apply to all spans)
        const udlResults = { v: Array(xShear.length).fill(0), m: Array(nNodes).fill(0), d: Array(nNodes).fill(0) };

        if (w_udl > 0) {
            // Single pass UDL on all spans
            const F_UDL = Array(nDOF).fill(0);
            let elemCounter = 0;
            for (const span of this.spans) {
                const le = span.length / nElemsPerSpan;
                const Fy = -w_udl * le / 2 * 1000;
                const Mom = -w_udl * le * le / 12 * 1000;

                for (let e = 0; e < nElemsPerSpan; e++) {
                    const nI = elemCounter * 2;
                    const nJ = (elemCounter + 1) * 2;
                    F_UDL[nI] += Fy;
                    F_UDL[nI + 1] += Mom;
                    F_UDL[nJ] += Fy;
                    F_UDL[nJ + 1] -= Mom;
                    elemCounter++;
                }
            }
            const U_UDL = this.solveSystem(K_Global, F_UDL, constrained);
            const forces = this.calculateForces(U_UDL, nTotalElems, elemLens, E, I, K_Global, F_UDL, constrained, nElemsPerSpan, numSpans);
            udlResults.v = forces.v;
            udlResults.m = forces.m;
            udlResults.d = forces.d;
        }

        // 6. Moving Truck Loop
        const truckLen = activeAxles.reduce((acc, a) => acc + a.spacing, 0);
        const totalLen = xNodes[xNodes.length - 1];
        const startPos = -truckLen;
        const endPos = totalLen + truckLen;

        // Helper to run a pass
        const runPass = (axles: Axle[]) => {
            for (let pos = startPos; pos <= endPos; pos += truckIncrement) {
                const F = Array(nDOF).fill(0);

                let axPos = pos;
                this.applyPointLoad(F, axPos, axles[0].load, elemLens, xNodes);

                for (let k = 0; k < axles.length - 1; k++) {
                    axPos -= axles[k].spacing;
                    this.applyPointLoad(F, axPos, axles[k + 1].load, elemLens, xNodes);
                }

                const U = this.solveSystem(K_Global, F, constrained);
                const { v, m, d, reactions } = this.calculateForces(U, nTotalElems, elemLens, E, I, K_Global, F, constrained, nElemsPerSpan, numSpans);

                // Update Envelopes
                // Shear (Stepped)
                for (let i = 0; i < v.length; i++) {
                    const val = v[i] + udlResults.v[i];
                    if (val > shearEnv[i].max) shearEnv[i].max = val;
                    if (val < shearEnv[i].min) shearEnv[i].min = val;
                }

                // Moment/Deflection (Nodal)
                for (let i = 0; i < m.length; i++) {
                    const valM = m[i] + udlResults.m[i];
                    if (valM > momEnv[i].max) momEnv[i].max = valM;
                    if (valM < momEnv[i].min) momEnv[i].min = valM;

                    const valD = d[i] + udlResults.d[i];
                    if (valD > defEnv[i].max) defEnv[i].max = valD;
                    if (valD < defEnv[i].min) defEnv[i].min = valD;
                }

                // Update Reaction Envelopes
                for (let i = 0; i < reactions.length; i++) {
                    if (reactions[i] > reactionEnv[i].max) reactionEnv[i].max = reactions[i];
                    if (reactions[i] < reactionEnv[i].min) reactionEnv[i].min = reactions[i];
                }
            }
        };

        // Forward Pass
        runPass(activeAxles);

        // Reverse Pass
        const revAxles = [...activeAxles].reverse().map((a, i, arr) => {
            const nextSpacing = i < arr.length - 1 ? arr[i + 1].spacing : 0;
            return { ...a, spacing: i < arr.length - 1 ? activeAxles[activeAxles.length - 2 - i].spacing : 0 };
        });

        // Fix spacings for reverse pass to match VBA logic exactly
        const vbaRevAxles = activeAxles.map((_, i) => activeAxles[activeAxles.length - 1 - i]);
        for (let i = 0; i < vbaRevAxles.length - 1; i++) {
            vbaRevAxles[i].spacing = activeAxles[activeAxles.length - 2 - i].spacing;
        }
        vbaRevAxles[vbaRevAxles.length - 1].spacing = 0;

        runPass(vbaRevAxles);

        return {
            shear: shearEnv,
            moment: momEnv,
            deflection: defEnv,
            xNodes,
            reactions: reactionEnv,
            supportPositions
        };
    }

    private applyPointLoad(F: number[], pos: number, mag: number, elemLens: number[], xNodes: number[]) {
        const totalLen = xNodes[xNodes.length - 1];
        if (pos < 0 || pos > totalLen) return;

        // Find element
        let elemIdx = -1;
        let localX = 0;

        // Quick search
        for (let i = 0; i < elemLens.length; i++) {
            if (pos <= xNodes[i + 1] + 0.000001) {
                elemIdx = i;
                localX = pos - xNodes[i];
                break;
            }
        }
        if (elemIdx === -1) elemIdx = elemLens.length - 1;

        const le = elemLens[elemIdx];
        const xi = localX / le;

        // Hermite Shape Functions
        const N1 = 1 - 3 * xi * xi + 2 * xi * xi * xi;
        const N2 = le * (xi - 2 * xi * xi + xi * xi * xi);
        const N3 = 3 * xi * xi - 2 * xi * xi * xi;
        const N4 = le * (-xi * xi + xi * xi * xi);

        const dofI = elemIdx * 2;
        const dofJ = (elemIdx + 1) * 2;

        F[dofI] -= mag * N1 * 1000;
        F[dofI + 1] -= mag * N2 * 1000;
        F[dofJ] -= mag * N3 * 1000;
        F[dofJ + 1] -= mag * N4 * 1000;
    }

    private calculateForces(
        U: number[],
        nElems: number,
        elemLens: number[],
        E: number,
        I: number,
        K_Global: number[][],
        F: number[],
        constrained: boolean[],
        nElemsPerSpan: number,
        numSpans: number
    ) {
        const v: number[] = [];
        const m: number[] = Array(nElems + 1).fill(0);
        const d: number[] = Array(nElems + 1).fill(0);
        const elemForces = [];

        // Element forces
        for (let e = 0; e < nElems; e++) {
            const le = elemLens[e];
            const coeff = (E * I) / Math.pow(le, 3);
            const u_loc = [U[e * 2], U[e * 2 + 1], U[(e + 1) * 2], U[(e + 1) * 2 + 1]];

            const k_loc = [
                [12, 6 * le, -12, 6 * le],
                [6 * le, 4 * le * le, -6 * le, 2 * le * le],
                [-12, -6 * le, 12, -6 * le],
                [6 * le, 2 * le * le, -6 * le, 4 * le * le]
            ];

            const f_loc = [0, 0, 0, 0];
            for (let r = 0; r < 4; r++) {
                for (let c = 0; c < 4; c++) {
                    f_loc[r] += k_loc[r][c] * coeff * u_loc[c];
                }
            }
            elemForces.push(f_loc);

            // Fill Shear (Stepped)
            v.push(f_loc[0] / 1000);      // Left edge
            v.push(-f_loc[2] / 1000);     // Right edge
        }

        // Nodal Averaging for Moment
        for (let n = 0; n <= nElems; n++) {
            d[n] = U[n * 2];

            let valM = 0;
            if (n === 0) {
                valM = elemForces[0][1];
            } else if (n === nElems) {
                valM = -elemForces[nElems - 1][3];
            } else {
                const M_left = -elemForces[n - 1][3];
                const M_right = elemForces[n][1];
                valM = (M_left + M_right) / 2;
            }
            m[n] = valM / 1000;
            m[n] = -m[n]; // Flip for convention
        }

        // Calculate Reactions at Supports (R = K*U - F at constrained DOFs)
        const nDOF = U.length;
        const reactions: number[] = [];

        // Get support node indices
        const supportNodes: number[] = [0]; // First support at node 0
        let nodeIdx = 0;
        for (let i = 0; i < numSpans; i++) {
            nodeIdx += nElemsPerSpan;
            supportNodes.push(nodeIdx);
        }

        // Calculate reaction at each support node
        for (const node of supportNodes) {
            const dof = node * 2; // Vertical DOF
            let reaction = 0;
            for (let j = 0; j < nDOF; j++) {
                reaction += K_Global[dof][j] * U[j];
            }
            reaction -= F[dof];
            reactions.push(-reaction / 1000); // Convert to kN, flip sign for upward positive
        }

        return { v, m, d, reactions };
    }
}

// --- COMPONENTS ---

// Helper to calculate nice ticks for chart axes
const calculateTicks = (min: number, max: number, targetCount: number) => {
    if (min === max) return [min];
    const span = max - min;
    const step = Math.pow(10, Math.floor(Math.log10(span / targetCount)));
    const err = targetCount / (span / step);

    let finalStep = step;
    if (err <= .15) finalStep *= 10;
    else if (err <= .35) finalStep *= 5;
    else if (err <= .75) finalStep *= 2;

    const start = Math.ceil(min / finalStep) * finalStep;
    const end = Math.floor(max / finalStep) * finalStep;

    const ticks = [];
    const decimals = Math.max(0, -Math.floor(Math.log10(finalStep)));

    for (let val = start; val <= end + (finalStep / 2); val += finalStep) {
        const cleanVal = parseFloat(val.toFixed(decimals));
        if (cleanVal >= min && cleanVal <= max) ticks.push(cleanVal);
    }
    return ticks;
};

const EnvelopeChart = ({
    data,
    dataKeyMax,
    dataKeyMin,
    title,
    unit,
    color,
    flipY = false
}: {
    data: any[],
    dataKeyMax: string,
    dataKeyMin: string,
    title: string,
    unit: string,
    color: string,
    flipY?: boolean
}) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const [width, setWidth] = useState(600);
    const height = 300;
    const padding = { top: 40, right: 30, bottom: 50, left: 70 };

    useEffect(() => {
        if (containerRef.current) setWidth(containerRef.current.clientWidth);
        const handleResize = () => containerRef.current && setWidth(containerRef.current.clientWidth);
        window.addEventListener('resize', handleResize);
        return () => window.removeEventListener('resize', handleResize);
    }, []);

    if (!data || data.length === 0) return <div className="h-[300px] flex items-center justify-center text-gray-400" > No Data </div>;

    const xVals = data.map(d => d.x);
    const maxVals = data.map(d => d[dataKeyMax]);
    const minVals = data.map(d => d[dataKeyMin]);
    const allY = [...maxVals, ...minVals];

    const xMin = Math.min(...xVals);
    const xMax = Math.max(...xVals);
    let yMin = Math.min(...allY);
    let yMax = Math.max(...allY);

    const yRange = yMax - yMin;
    if (yRange === 0) { yMax += 1; yMin -= 1; }
    else { yMax += yRange * 0.1; yMin -= yRange * 0.1; }

    const xTicks = calculateTicks(xMin, xMax, 8);
    const yTicks = calculateTicks(yMin, yMax, 6);

    const xScale = (val: number) => padding.left + ((val - xMin) / (xMax - xMin)) * (width - padding.left - padding.right);
    const yScale = (val: number) => height - padding.bottom - ((val - yMin) / (yMax - yMin)) * (height - padding.top - padding.bottom);

    const createPath = (vals: number[]) => {
        return vals.map((y, i) => `${i === 0 ? 'M' : 'L'} ${xScale(data[i].x)} ${yScale(y)}`).join(' ');
    };

    const pathMax = createPath(maxVals);
    const pathMin = createPath(minVals);

    const pathFill = `${pathMax} L ${xScale(data[data.length - 1].x)} ${yScale(minVals[minVals.length - 1])} ` +
        minVals.slice().reverse().map((y, i) => `L ${xScale(data[data.length - 1 - i].x)} ${yScale(y)}`).join(' ') + " Z";

    const zeroY = yScale(0);

    return (
        <div ref={containerRef} className="w-full bg-white rounded-lg shadow-sm border border-gray-200 p-4 mb-6" >
            <h3 className="text-lg font-semibold text-gray-800 mb-2" > {title} </h3>
            < svg width={width} height={height} className="overflow-visible" >

                {/* X-Grid & Labels */}
                {
                    xTicks.map(tick => {
                        const xPos = xScale(tick);
                        return (
                            <g key={`x-${tick}`
                            }>
                                <line x1={xPos} y1={padding.top} x2={xPos} y2={height - padding.bottom} stroke="#e5e7eb" strokeWidth="1" />
                                <text x={xPos} y={height - padding.bottom + 15
                                } textAnchor="middle" fontSize="10" fill="#6b7280" > {tick} </text>
                            </g>
                        );
                    })}

                {/* Y-Grid & Labels */}
                {
                    yTicks.map(tick => {
                        const yPos = yScale(tick);
                        return (
                            <g key={`y-${tick}`
                            }>
                                <line x1={padding.left} y1={yPos} x2={width - padding.right} y2={yPos} stroke="#e5e7eb" strokeWidth="1" />
                                <text x={padding.left - 8} y={yPos + 3} textAnchor="end" fontSize="10" fill="#6b7280" > {tick} </text>
                            </g>
                        );
                    })}

                <line x1={padding.left} y1={padding.top} x2={padding.left} y2={height - padding.bottom} stroke="#374151" strokeWidth="1" />
                <line x1={padding.left} y1={height - padding.bottom} x2={width - padding.right} y2={height - padding.bottom} stroke="#374151" strokeWidth="1" />

                {zeroY > padding.top && zeroY < height - padding.bottom && (
                    <line x1={padding.left} y1={zeroY} x2={width - padding.right} y2={zeroY} stroke="#9ca3af" strokeWidth="1.5" strokeDasharray="4 4" />
                )}

                <text
                    x={padding.left + (width - padding.left - padding.right) / 2}
                    y={height - 10}
                    textAnchor="middle"
                    fontSize="12"
                    fontWeight="500"
                    fill="#374151"
                >
                    Length(m)
                </text>

                < text
                    x={15}
                    y={padding.top + (height - padding.top - padding.bottom) / 2}
                    textAnchor="middle"
                    fontSize="12"
                    fontWeight="500"
                    fill="#374151"
                    className="transform -rotate-90 origin-center"
                    style={{ transformBox: 'fill-box' }}
                >
                    {unit}
                </text>

                < path d={pathFill} fill={color} fillOpacity="0.1" />
                <path d={pathMax} fill="none" stroke={color} strokeWidth="2" />
                <path d={pathMin} fill="none" stroke="#ef4444" strokeWidth="2" strokeDasharray="4 2" />

            </svg>
            < div className="flex justify-center gap-6 mt-2 text-sm" >
                <div className="flex items-center" > <div className="w-4 h-0.5 bg-[color:var(--color)] mr-2" style={{ backgroundColor: color }}> </div> Max Envelope</div >
                <div className="flex items-center" > <div className="w-4 h-0.5 bg-red-500 mr-2 border-dashed border-t-2 border-red-500" > </div> Min Envelope</div >
            </div>
        </div>
    );
};

// Beam Schematic for Configuration View (shows beam with supports)
const BeamSchematic = ({ spans }: { spans: Span[] }) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const [width, setWidth] = useState(600);
    const height = 120;
    const padding = { top: 30, right: 40, bottom: 40, left: 40 };

    useEffect(() => {
        if (containerRef.current) setWidth(containerRef.current.clientWidth);
        const handleResize = () => containerRef.current && setWidth(containerRef.current.clientWidth);
        window.addEventListener('resize', handleResize);
        return () => window.removeEventListener('resize', handleResize);
    }, []);

    const totalLen = spans.reduce((a, b) => a + b.length, 0);
    if (totalLen === 0) return null;

    const beamY = height / 2;
    const xScale = (val: number) => padding.left + (val / totalLen) * (width - padding.left - padding.right);

    // Calculate support positions
    const supportPositions: number[] = [0];
    let cumLen = 0;
    for (const span of spans) {
        cumLen += span.length;
        supportPositions.push(cumLen);
    }

    // Triangle support symbol
    const triangleSize = 12;
    const renderSupport = (x: number, idx: number) => {
        const xPos = xScale(x);
        return (
            <g key={`support-${idx}`}>
                {/* Triangle */}
                <polygon
                    points={`${xPos},${beamY + 4} ${xPos - triangleSize},${beamY + triangleSize + 8} ${xPos + triangleSize},${beamY + triangleSize + 8}`}
                    fill="#374151"
                    stroke="#1f2937"
                    strokeWidth="1"
                />
                {/* Ground line */}
                <line
                    x1={xPos - triangleSize - 4}
                    y1={beamY + triangleSize + 10}
                    x2={xPos + triangleSize + 4}
                    y2={beamY + triangleSize + 10}
                    stroke="#1f2937"
                    strokeWidth="2"
                />
                {/* Support label */}
                <text
                    x={xPos}
                    y={beamY + triangleSize + 25}
                    textAnchor="middle"
                    fontSize="10"
                    fill="#6b7280"
                >
                    {idx === 0 ? 'A' : String.fromCharCode(65 + idx)}
                </text>
            </g>
        );
    };

    return (
        <div ref={containerRef} className="w-full bg-white rounded-lg shadow-sm border border-gray-200 p-4 mb-4">
            <h3 className="text-sm font-semibold text-gray-700 mb-2">Beam Configuration</h3>
            <svg width={width} height={height} className="overflow-visible">
                {/* Beam line */}
                <line
                    x1={xScale(0)}
                    y1={beamY}
                    x2={xScale(totalLen)}
                    y2={beamY}
                    stroke="#2563eb"
                    strokeWidth="6"
                    strokeLinecap="round"
                />

                {/* Span labels and dimension lines */}
                {spans.map((span, idx) => {
                    let startX = 0;
                    for (let i = 0; i < idx; i++) startX += spans[i].length;
                    const endX = startX + span.length;
                    const midX = (startX + endX) / 2;

                    return (
                        <g key={`span-${idx}`}>
                            {/* Dimension line */}
                            <line
                                x1={xScale(startX) + 2}
                                y1={beamY - 20}
                                x2={xScale(endX) - 2}
                                y2={beamY - 20}
                                stroke="#9ca3af"
                                strokeWidth="1"
                                markerStart="url(#arrowLeft)"
                                markerEnd="url(#arrowRight)"
                            />
                            {/* Span length label */}
                            <text
                                x={xScale(midX)}
                                y={beamY - 25}
                                textAnchor="middle"
                                fontSize="11"
                                fill="#374151"
                                fontWeight="500"
                            >
                                {span.length}m
                            </text>
                        </g>
                    );
                })}

                {/* Arrow markers definition */}
                <defs>
                    <marker id="arrowLeft" markerWidth="6" markerHeight="6" refX="0" refY="3" orient="auto">
                        <path d="M6,0 L0,3 L6,6" fill="none" stroke="#9ca3af" strokeWidth="1" />
                    </marker>
                    <marker id="arrowRight" markerWidth="6" markerHeight="6" refX="6" refY="3" orient="auto">
                        <path d="M0,0 L6,3 L0,6" fill="none" stroke="#9ca3af" strokeWidth="1" />
                    </marker>
                </defs>

                {/* Supports */}
                {supportPositions.map((pos, idx) => renderSupport(pos, idx))}
            </svg>
        </div>
    );
};

// Beam Reaction Diagram for Results View (shows beam with max reaction values)
const BeamReactionDiagram = ({
    spans,
    reactions,
    supportPositions
}: {
    spans: Span[],
    reactions: ReactionEnvelope[],
    supportPositions: number[]
}) => {
    const containerRef = useRef<HTMLDivElement>(null);
    const [width, setWidth] = useState(600);
    const height = 200;
    const padding = { top: 30, right: 60, bottom: 60, left: 60 };

    useEffect(() => {
        if (containerRef.current) setWidth(containerRef.current.clientWidth);
        const handleResize = () => containerRef.current && setWidth(containerRef.current.clientWidth);
        window.addEventListener('resize', handleResize);
        return () => window.removeEventListener('resize', handleResize);
    }, []);

    const totalLen = spans.reduce((a, b) => a + b.length, 0);
    if (totalLen === 0 || !reactions || reactions.length === 0) return null;

    const beamY = 55;
    const xScale = (val: number) => padding.left + (val / totalLen) * (width - padding.left - padding.right);

    const triangleSize = 10;

    const renderSupportWithReaction = (x: number, idx: number, reaction: ReactionEnvelope) => {
        const xPos = xScale(x);
        const maxR = Math.abs(reaction.max);
        const minR = Math.abs(reaction.min);
        const displayR = Math.max(maxR, minR);

        return (
            <g key={`support-${idx}`}>
                {/* Triangle support */}
                <polygon
                    points={`${xPos},${beamY + 4} ${xPos - triangleSize},${beamY + triangleSize + 6} ${xPos + triangleSize},${beamY + triangleSize + 6}`}
                    fill="#374151"
                    stroke="#1f2937"
                    strokeWidth="1"
                />

                {/* Ground line under triangle */}
                <line
                    x1={xPos - triangleSize - 3}
                    y1={beamY + triangleSize + 8}
                    x2={xPos + triangleSize + 3}
                    y2={beamY + triangleSize + 8}
                    stroke="#1f2937"
                    strokeWidth="2"
                />

                {/* Support label (A, B, C, D) - positioned to the side */}
                <text
                    x={xPos}
                    y={beamY + triangleSize + 25}
                    textAnchor="middle"
                    fontSize="11"
                    fill="#374151"
                    fontWeight="600"
                >
                    {String.fromCharCode(65 + idx)}
                </text>

                {/* Reaction arrow (upward) - line and arrowhead */}
                <line
                    x1={xPos}
                    y1={beamY + triangleSize + 55}
                    x2={xPos}
                    y2={beamY + triangleSize + 35}
                    stroke="#dc2626"
                    strokeWidth="2.5"
                />
                {/* Arrowhead pointing UP */}
                <polygon
                    points={`${xPos - 5},${beamY + triangleSize + 35} ${xPos},${beamY + triangleSize + 27} ${xPos + 5},${beamY + triangleSize + 35}`}
                    fill="#dc2626"
                />

                {/* Reaction value - at bottom */}
                <text
                    x={xPos}
                    y={beamY + triangleSize + 72}
                    textAnchor="middle"
                    fontSize="11"
                    fill="#dc2626"
                    fontWeight="600"
                >
                    {displayR.toFixed(1)} kN
                </text>
            </g>
        );
    };

    return (
        <div ref={containerRef} className="w-full bg-white rounded-lg shadow-sm border border-gray-200 p-4 mb-6">
            <h3 className="text-lg font-semibold text-gray-800 mb-2">Maximum Support Reactions</h3>
            <svg width={width} height={height} className="overflow-visible">
                {/* Reaction arrow marker - pointing upward */}
                <defs>
                    <marker id="reactionArrow" markerWidth="10" markerHeight="10" refX="5" refY="10" orient="auto">
                        <path d="M0,10 L5,0 L10,10 Z" fill="#dc2626" />
                    </marker>
                </defs>

                {/* Beam line */}
                <line
                    x1={xScale(0)}
                    y1={beamY}
                    x2={xScale(totalLen)}
                    y2={beamY}
                    stroke="#2563eb"
                    strokeWidth="6"
                    strokeLinecap="round"
                />

                {/* Span labels - above beam */}
                {spans.map((span, idx) => {
                    let startX = 0;
                    for (let i = 0; i < idx; i++) startX += spans[i].length;
                    const endX = startX + span.length;
                    const midX = (startX + endX) / 2;

                    return (
                        <g key={`span-label-${idx}`}>
                            {/* Dimension line */}
                            <line
                                x1={xScale(startX) + 5}
                                y1={beamY - 18}
                                x2={xScale(endX) - 5}
                                y2={beamY - 18}
                                stroke="#9ca3af"
                                strokeWidth="1"
                            />
                            {/* End ticks */}
                            <line x1={xScale(startX) + 5} y1={beamY - 14} x2={xScale(startX) + 5} y2={beamY - 22} stroke="#9ca3af" strokeWidth="1" />
                            <line x1={xScale(endX) - 5} y1={beamY - 14} x2={xScale(endX) - 5} y2={beamY - 22} stroke="#9ca3af" strokeWidth="1" />
                            {/* Label */}
                            <text
                                x={xScale(midX)}
                                y={beamY - 25}
                                textAnchor="middle"
                                fontSize="10"
                                fill="#6b7280"
                            >
                                {span.length}m
                            </text>
                        </g>
                    );
                })}

                {/* Supports with reactions */}
                {supportPositions.map((pos, idx) =>
                    reactions[idx] && renderSupportWithReaction(pos, idx, reactions[idx])
                )}
            </svg>
            <div className="flex justify-center gap-4 mt-1 text-xs text-gray-500">
                <span>↑ Maximum reaction forces shown</span>
            </div>
        </div>
    );
};

export default function BeamAnalysisApp() {
    const [spans, setSpans] = useState<Span[]>(DEFAULT_SPANS);
    const [axles, setAxles] = useState<Axle[]>(DEFAULT_AXLES);
    const [config, setConfig] = useState<AnalysisConfig>(DEFAULT_CONFIG);
    const [results, setResults] = useState<AnalysisResults | null>(null);
    const [isAnalyzing, setIsAnalyzing] = useState(false);
    const [activeTab, setActiveTab] = useState<'config' | 'results'>('config');

    // Inject SheetJS script for Excel Export
    useEffect(() => {
        const script = document.createElement('script');
        script.src = "https://cdn.sheetjs.com/xlsx-0.20.1/package/dist/xlsx.full.min.js";
        script.async = true;
        document.body.appendChild(script);
        return () => {
            document.body.removeChild(script);
        }
    }, []);

    const addSpan = () => {
        const newId = `s${Date.now()}`;
        setSpans([...spans, { id: newId, length: 20 }]);
    };

    const removeSpan = (id: string) => {
        if (spans.length <= 1) return;
        setSpans(spans.filter(s => s.id !== id));
    };

    const updateSpan = (id: string, val: number) => {
        setSpans(spans.map(s => s.id === id ? { ...s, length: val } : s));
    };

    const runAnalysis = async () => {
        setIsAnalyzing(true);
        setTimeout(() => {
            try {
                const solver = new BeamFEM(spans, axles, config);
                const res = solver.runAnalysis();
                setResults(res);
                setActiveTab('results');
            } catch (e) {
                console.error(e);
                alert("Analysis failed. Check inputs.");
            }
            setIsAnalyzing(false);
        }, 100);
    };

    const downloadExcel = () => {
        if (!results) return;

        // @ts-ignore
        if (typeof window === 'undefined' || !window.XLSX) {
            alert("Excel export library is loading. Please try again in a few seconds.");
            return;
        }

        // @ts-ignore
        const wb = window.XLSX.utils.book_new();

        const formatData = (data: EnvelopePoint[]) => data.map(d => ({
            "Position (m)": parseFloat(d.x.toFixed(3)),
            "Max": parseFloat(d.max.toFixed(3)),
            "Min": parseFloat(d.min.toFixed(3))
        }));

        // Shear Sheet
        // @ts-ignore
        const wsShear = window.XLSX.utils.json_to_sheet(formatData(results.shear));
        // @ts-ignore
        window.XLSX.utils.book_append_sheet(wb, wsShear, "Shear Force");

        // Moment Sheet
        // @ts-ignore
        const wsMoment = window.XLSX.utils.json_to_sheet(formatData(results.moment));
        // @ts-ignore
        window.XLSX.utils.book_append_sheet(wb, wsMoment, "Bending Moment");

        // Deflection Sheet
        // @ts-ignore
        const wsDef = window.XLSX.utils.json_to_sheet(formatData(results.deflection));
        // @ts-ignore
        window.XLSX.utils.book_append_sheet(wb, wsDef, "Deflection");

        // Write file
        // @ts-ignore
        window.XLSX.writeFile(wb, "beam_analysis_results.xlsx");
    };

    return (
        <div className="min-h-screen bg-gray-50 text-slate-800 font-sans" >
            {/* Header */}
            < header className="bg-blue-700 text-white p-4 shadow-md sticky top-0 z-10" >
                <div className="max-w-5xl mx-auto flex justify-between items-center" >
                    <h1 className="text-xl font-bold flex items-center gap-2" >
                        <span className="bg-white text-blue-700 p-1 rounded font-black text-xs" > FEM </span>
                        Beam Analysis < span className="text-blue-200 font-normal text-sm hidden sm:inline" >| CL - 625 Truck Moving Load </span>
                    </h1>
                    < button
                        onClick={runAnalysis}
                        disabled={isAnalyzing}
                        className={`flex items-center gap-2 px-4 py-2 rounded font-medium transition-colors ${isAnalyzing ? 'bg-blue-800 cursor-wait' : 'bg-white text-blue-700 hover:bg-blue-50'}`
                        }
                    >
                        {isAnalyzing ? <RotateCcw className="animate-spin w-4 h-4" /> : <Play className="w-4 h-4" />}
                        {isAnalyzing ? 'Calculating...' : 'Run Analysis'}
                    </button>
                </div>
            </header>

            < main className="max-w-5xl mx-auto p-4 md:p-6" >

                {/* Tabs */}
                < div className="flex gap-4 border-b border-gray-200 mb-6" >
                    <button
                        onClick={() => setActiveTab('config')}
                        className={`pb-2 px-1 font-medium text-sm transition-colors ${activeTab === 'config' ? 'text-blue-600 border-b-2 border-blue-600' : 'text-gray-500 hover:text-gray-700'}`}
                    >
                        Configuration
                    </button>
                    < button
                        onClick={() => setActiveTab('results')}
                        disabled={!results}
                        className={`pb-2 px-1 font-medium text-sm transition-colors ${activeTab === 'results' ? 'text-blue-600 border-b-2 border-blue-600' : 'text-gray-500 hover:text-gray-700 disabled:opacity-50'}`}
                    >
                        Results
                    </button>
                </div>

                {
                    activeTab === 'config' && (
                        <>
                            <BeamSchematic spans={spans} />
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-6" >
                                {/* Spans Card */}
                                < div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200" >
                                    <div className="flex justify-between items-center mb-4" >
                                        <h2 className="text-lg font-semibold flex items-center gap-2" >
                                            <Settings className="w-5 h-5 text-gray-500" /> Geometry
                                        </h2>
                                        < button onClick={addSpan} className="text-sm bg-blue-50 text-blue-600 px-3 py-1 rounded hover:bg-blue-100 flex items-center gap-1" >
                                            <Plus className="w-3 h-3" /> Add Span
                                        </button>
                                    </div>

                                    < div className="space-y-3" >
                                        {
                                            spans.map((span, idx) => (
                                                <div key={span.id} className="flex items-center gap-3 p-3 bg-gray-50 rounded border border-gray-100" >
                                                    <span className="text-sm font-bold text-gray-400 w-8" >#{idx + 1} </span>
                                                    < div className="flex-1" >
                                                        <label className="text-xs text-gray-500 block" > Length(m) </label>
                                                        < input
                                                            type="number"
                                                            value={span.length}
                                                            onChange={(e) => updateSpan(span.id, parseFloat(e.target.value) || 0)
                                                            }
                                                            className="w-full bg-white border border-gray-300 rounded px-2 py-1 text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                                                        />
                                                    </div>
                                                    < button onClick={() => removeSpan(span.id)} className="text-gray-400 hover:text-red-500 p-2" >
                                                        <Trash2 className="w-4 h-4" />
                                                    </button>
                                                </div>
                                            ))}
                                    </div>
                                    < div className="mt-4 pt-4 border-t border-gray-100 text-sm text-gray-500 flex justify-between" >
                                        <span>Total Length: </span>
                                        < span className="font-mono font-bold text-gray-800" > {spans.reduce((a, b) => a + b.length, 0).toFixed(2)} m </span>
                                    </div>
                                </div>

                                {/* Config Card */}
                                <div className="space-y-6" >
                                    <div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200" >
                                        <h2 className="text-lg font-semibold mb-4" > Analysis Settings </h2>
                                        < div className="space-y-4" >
                                            <div>
                                                <label className="block text-sm font-medium text-gray-700 mb-1" > Load Case </label>
                                                < select
                                                    value={config.loadCase}
                                                    onChange={(e) => setConfig({ ...config, loadCase: e.target.value as 'truck' | 'lane' })}
                                                    className="w-full bg-white border border-gray-300 rounded px-3 py-2 text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                                                >
                                                    <option value="truck" > CL - 625 Truck Only(Standard) </option>
                                                    < option value="lane" > CL - 625 Lane Load(80 % Truck + 9 kN / m) </option>
                                                </select>
                                            </div>

                                            < div className="grid grid-cols-2 gap-4" >
                                                <div>
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Elastic Modulus(Pa) </label>
                                                    < input
                                                        type="number"
                                                        value={config.E}
                                                        onChange={(e) => setConfig({ ...config, E: parseFloat(e.target.value) })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                                < div >
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Inertia(m⁴) </label>
                                                    < input
                                                        type="number"
                                                        value={config.I}
                                                        onChange={(e) => setConfig({ ...config, I: parseFloat(e.target.value) })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                            </div>

                                            < div className="bg-blue-50 p-3 rounded text-sm text-blue-800 flex gap-2 items-start" >
                                                <AlertCircle className="w-4 h-4 mt-0.5 shrink-0" />
                                                <p>Mesh refinement is set to {config.nElemsPerSpan} elements per span.Truck moves in {config.truckIncrement}m increments.</p>
                                            </div>
                                        </div>
                                    </div>

                                    < div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200" >
                                        <h2 className="text-lg font-semibold mb-2" > Truck Configuration </h2>
                                        < p className="text-sm text-gray-500 mb-4" > Standard CL - 625 Axle Loads(kN) and Spacings(m) </p>
                                        < div className="flex flex-wrap gap-2" >
                                            {
                                                axles.map((axle, i) => (
                                                    <div key={axle.id} className="bg-gray-100 rounded p-2 text-center min-w-[60px]" >
                                                        <div className="text-xs font-bold text-gray-500" > Axle {i + 1} </div>
                                                        < div className="font-mono text-sm text-blue-600 font-bold" > {axle.load} </div>
                                                        {
                                                            i < axles.length - 1 && (
                                                                <div className="text-[10px] text-gray-400 mt-1 border-t border-gray-300 pt-1" >
                                                                    ↓ {axle.spacing} m
                                                                </div>
                                                            )
                                                        }
                                                    </div>
                                                ))}
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </>
                    )}

                {
                    activeTab === 'results' && results && (
                        <div className="animate-in fade-in slide-in-from-bottom-4 duration-500" >
                            <BeamReactionDiagram
                                spans={spans}
                                reactions={results.reactions}
                                supportPositions={results.supportPositions}
                            />
                            <EnvelopeChart
                                title="Shear Force Envelope"
                                data={results.shear}
                                dataKeyMax="max"
                                dataKeyMin="min"
                                unit="Shear (kN)"
                                color="#2563eb"
                            />

                            <EnvelopeChart
                                title="Bending Moment Envelope"
                                data={results.moment}
                                dataKeyMax="max"
                                dataKeyMin="min"
                                unit="Moment (kNm)"
                                color="#059669"
                                flipY={true}
                            />

                            <EnvelopeChart
                                title="Deflection Envelope"
                                data={results.deflection}
                                dataKeyMax="max"
                                dataKeyMin="min"
                                unit="Deflection (m)"
                                color="#9333ea"
                            />

                            <div className="bg-white p-6 rounded-lg shadow-sm border border-gray-200 mt-6" >
                                <h3 className="font-semibold mb-4" > Export Data </h3>
                                < div className="text-sm text-gray-600 mb-4" >
                                    Download the analysis results as an Excel(.xlsx) file with separate sheets for Shear, Moment, and Deflection.
                                </div>
                                < button
                                    onClick={downloadExcel}
                                    className="bg-green-700 text-white px-4 py-2 rounded text-sm hover:bg-green-800 flex items-center gap-2"
                                >
                                    <Download className="w-4 h-4" />
                                    Download Excel(XLSX)
                                </button>
                            </div>
                        </div>
                    )
                }
            </main>
        </div>
    );
}