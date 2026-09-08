import { useState, useEffect, useRef, useMemo } from 'react';
import { Play, RotateCcw, Plus, Trash2, Settings, AlertCircle, Download } from 'lucide-react';

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
    loadCase: 'truck' | 'lane' | 'envelope';
    dlaOverride?: number | null;
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
    govPos?: number; // truck lead position producing max reaction
};

type AnalysisResults = {
    shear: EnvelopePoint[];
    moment: EnvelopePoint[];
    deflection: EnvelopePoint[];
    xNodes: number[];
    reactions: ReactionEnvelope[];
    supportPositions: number[];
    dlaUsed?: number;
    incrementUsed?: number;
    baseIncrement?: number;
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

const MAX_PATTERN_SPANS = 12; // 2^n UDL patterns becomes prohibitive beyond this
const MIN_SWEEP_STEP = 0.02; // m - floor to avoid runaway solve counts on tiny spans
const MAX_SWEEP_STEPS = 6000; // per pass - cap to keep UI responsive on very long bridges

// Adaptive moving-load step. A fixed 0.25m step can jump over entire elements on
// short trestle spans (e.g. 4m/40 = 0.10m elements) and miss peaks between supports.
// Effective step = min(base, minSpan/40, minElemLen, minAxleSpacing/8), floored and
// capped so total uniform steps stay <= MAX_SWEEP_STEPS. Support alignments are
// always added explicitly, so support peaks are exact regardless of step.
export function computeEffectiveIncrement(
    spans: Span[],
    axles: Axle[],
    baseIncrement: number,
    nElemsPerSpan: number
): { effective: number; reason: string; wasReduced: boolean } {
    const lens = spans.map(s => s.length).filter(l => Number.isFinite(l) && l > 0);
    if (lens.length === 0 || !Number.isFinite(baseIncrement) || baseIncrement <= 0) {
        return { effective: baseIncrement, reason: '', wasReduced: false };
    }
    const minSpan = Math.min(...lens);
    const totalLen = lens.reduce((a, b) => a + b, 0);
    const minElem = minSpan / Math.max(2, nElemsPerSpan);
    const gaps = axles.slice(0, -1).map(a => a.spacing).filter(s => Number.isFinite(s) && s > 0);
    const minGap = gaps.length > 0 ? Math.min(...gaps) : Infinity;
    const truckLen = axles.reduce((acc, a) => acc + (Number.isFinite(a.spacing) ? a.spacing : 0), 0);

    let eff = baseIncrement;
    let control = 'base setting';
    const consider = (cand: number, label: string) => {
        if (Number.isFinite(cand) && cand > 0 && cand < eff) {
            eff = cand;
            control = label;
        }
    };
    consider(minSpan / 40, `short span control (minSpan ${minSpan.toFixed(2)}m / 40)`);
    consider(minElem, `element control (minElem ${minElem.toFixed(3)}m)`);
    if (Number.isFinite(minGap)) consider(minGap / 8, `axle-spacing control (minGap ${minGap.toFixed(2)}m / 8)`);

    if (eff < MIN_SWEEP_STEP) {
        eff = MIN_SWEEP_STEP;
        control = `minimum step floor (${MIN_SWEEP_STEP.toFixed(2)}m)`;
    }
    // Cap total uniform steps for very long bridges / tiny steps
    const sweepLen = totalLen + 2 * truckLen;
    if (sweepLen / eff > MAX_SWEEP_STEPS) {
        eff = sweepLen / MAX_SWEEP_STEPS;
        control = `step-count cap (${MAX_SWEEP_STEPS} steps over ${sweepLen.toFixed(1)}m sweep)`;
    }
    // Round up to 5mm to avoid pathological 0.09333... steps; rounding up never refines beyond need
    eff = Math.ceil(eff * 200) / 200;
    if (eff > baseIncrement) eff = baseIncrement;
    const wasReduced = eff < baseIncrement - 1e-9;
    return { effective: eff, reason: wasReduced ? control : '', wasReduced };
}

const DEFAULT_CONFIG: AnalysisConfig = {
    E: 200000000000,
    I: 0.005,
    nElemsPerSpan: 40,
    truckIncrement: 0.25,
    loadCase: 'truck',
    dlaOverride: null,
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

    private static maxAxlesOnLength(spanLen: number, axles: Axle[]): number {
        let maxCount = 1;
        for (let i = 0; i < axles.length; i++) {
            let count = 1;
            let cumulativeDist = 0;
            for (let j = i; j < axles.length - 1; j++) {
                cumulativeDist += axles[j].spacing;
                if (cumulativeDist <= spanLen + 0.000001) {
                    count++;
                } else {
                    break;
                }
            }
            if (count > maxCount) maxCount = count;
        }
        return maxCount;
    }

    private static dlaForAxleCount(n: number): number {
        if (n <= 1) return 0.40;
        if (n === 2) return 0.30;
        return 0.25;
    }

    private calculateAutoDla(): number {
        // Conservative: evaluate DLA per span and take the maximum.
        // Using only maxSpan would under-apply DLA to short spans.
        let worst = 0.25;
        for (const span of this.spans) {
            const n = BeamFEM.maxAxlesOnLength(span.length, this.axles);
            const d = BeamFEM.dlaForAxleCount(n);
            if (d > worst) worst = d;
        }
        return worst;
    }

    public runAnalysis(): AnalysisResults {
        const { E, I, nElemsPerSpan, truckIncrement, loadCase } = this.config;
        // --- Input validation (mirrors VBA MainAnalysis guards) ---
        if (!this.spans || this.spans.length < 1) throw new Error('At least one span is required.');
        for (let i = 0; i < this.spans.length; i++) {
            if (!Number.isFinite(this.spans[i].length) || this.spans[i].length <= 0)
                throw new Error(`Span ${i + 1} length must be greater than zero.`);
        }
        if (!Number.isFinite(E) || E <= 0) throw new Error('Elastic modulus E must be greater than zero.');
        if (!Number.isFinite(I) || I <= 0) throw new Error('Moment of inertia I must be greater than zero.');
        if (!Number.isFinite(nElemsPerSpan) || nElemsPerSpan < 2) throw new Error('Elements per span must be >= 2.');
        if (!Number.isFinite(truckIncrement) || truckIncrement <= 0) throw new Error('Truck increment must be > 0.');
        if (!this.axles || this.axles.length < 1) throw new Error('At least one axle is required.');
        for (let i = 0; i < this.axles.length; i++) {
            if (!Number.isFinite(this.axles[i].load) || this.axles[i].load < 0)
                throw new Error(`Axle ${i + 1} load must be >= 0.`);
            if (i < this.axles.length - 1 && (!Number.isFinite(this.axles[i].spacing) || this.axles[i].spacing < 0))
                throw new Error(`Axle ${i + 1} spacing must be >= 0.`);
        }
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

        // 4. Envelope factory (shear: 2 pts/elem stepped; moment/deflection: nodal)
        const initEnvelope = (len: number) => Array(len).fill(0).map(() => ({ max: -Infinity, min: Infinity }));

        // Reaction envelopes (one per support) - filled by runSingleCase()
        const autoDla = this.calculateAutoDla();
        const dla = (this.config.dlaOverride !== undefined && this.config.dlaOverride !== null)
            ? this.config.dlaOverride
            : autoDla;
        if (!Number.isFinite(dla) || dla < 0) throw new Error('DLA must be >= 0.');

        const baseAxles: Axle[] = this.axles.map(a => ({ ...a })); // deep copy, never mutate this.axles
        const spanLengths = this.spans.map(s => s.length);
        // Adaptive sweep step: fixed 0.25m skips elements on short trestle spans.
        const { effective: step } = computeEffectiveIncrement(this.spans, this.axles, truckIncrement, nElemsPerSpan);

        type UdlEnvs = {
            Vmax: number[]; Vmin: number[];
            Mmax: number[]; Mmin: number[];
            Dmax: number[]; Dmin: number[];
            Rmax: number[]; Rmin: number[];
        };

        // Patterned lane-UDL envelopes (2^n span patterns, VBA CalculateUDLEnvelopes parity).
        // All-spans-loaded alone misses alternate-span governing maxima, so enumerate patterns.
        const calculateUdlEnvelopes = (w_udl: number): UdlEnvs => {
            const nShear = nTotalElems * 2;
            const Vmax = Array(nShear).fill(-Infinity);
            const Vmin = Array(nShear).fill(Infinity);
            const Mmax = Array(nNodes).fill(-Infinity);
            const Mmin = Array(nNodes).fill(Infinity);
            const Dmax = Array(nNodes).fill(-Infinity);
            const Dmin = Array(nNodes).fill(Infinity);
            const Rmax = Array(supportPositions.length).fill(-Infinity);
            const Rmin = Array(supportPositions.length).fill(Infinity);
            if (!(w_udl > 0)) {
                return {
                    Vmax: Array(nShear).fill(0), Vmin: Array(nShear).fill(0),
                    Mmax: Array(nNodes).fill(0), Mmin: Array(nNodes).fill(0),
                    Dmax: Array(nNodes).fill(0), Dmin: Array(nNodes).fill(0),
                    Rmax: Array(supportPositions.length).fill(0), Rmin: Array(supportPositions.length).fill(0),
                };
            }
            if (numSpans > MAX_PATTERN_SPANS)
                throw new Error(`Lane load patterning evaluates 2^n patterns and is limited to ${MAX_PATTERN_SPANS} spans.`);
            const nPatterns = 1 << numSpans;
            for (let p = 0; p < nPatterns; p++) {
                const F_UDL = Array(nDOF).fill(0);
                const elemLoads_UDL = Array.from({ length: nTotalElems }, () => [0, 0, 0, 0]);
                let elemCounter = 0;
                for (let s = 0; s < numSpans; s++) {
                    const isLoaded = ((p >> s) & 1) === 1;
                    const le = spanLengths[s] / nElemsPerSpan;
                    if (isLoaded) {
                        const Fy = -w_udl * le / 2 * 1000;
                        const Mom = -w_udl * le * le / 12 * 1000;
                        for (let e = 0; e < nElemsPerSpan; e++) {
                            const nI = elemCounter * 2;
                            const nJ = (elemCounter + 1) * 2;
                            F_UDL[nI] += Fy;
                            F_UDL[nI + 1] += Mom;
                            F_UDL[nJ] += Fy;
                            F_UDL[nJ + 1] -= Mom;
                            elemLoads_UDL[elemCounter][0] = Fy;
                            elemLoads_UDL[elemCounter][1] = Mom;
                            elemLoads_UDL[elemCounter][2] = Fy;
                            elemLoads_UDL[elemCounter][3] = -Mom;
                            elemCounter++;
                        }
                    } else {
                        elemCounter += nElemsPerSpan;
                    }
                }
                const U_UDL = solveFactored(F_UDL);
                const forces = this.calculateForces(U_UDL, nTotalElems, elemLens, E, I, K_Global, F_UDL, constrained, nElemsPerSpan, numSpans, elemLoads_UDL);
                for (let i = 0; i < forces.v.length; i++) {
                    if (forces.v[i] > Vmax[i]) Vmax[i] = forces.v[i];
                    if (forces.v[i] < Vmin[i]) Vmin[i] = forces.v[i];
                }
                for (let i = 0; i < forces.m.length; i++) {
                    if (forces.m[i] > Mmax[i]) Mmax[i] = forces.m[i];
                    if (forces.m[i] < Mmin[i]) Mmin[i] = forces.m[i];
                    if (forces.d[i] > Dmax[i]) Dmax[i] = forces.d[i];
                    if (forces.d[i] < Dmin[i]) Dmin[i] = forces.d[i];
                }
                for (let s = 0; s < forces.reactions.length; s++) {
                    if (forces.reactions[s] > Rmax[s]) Rmax[s] = forces.reactions[s];
                    if (forces.reactions[s] < Rmin[s]) Rmin[s] = forces.reactions[s];
                }
            }
            return { Vmax, Vmin, Mmax, Mmin, Dmax, Dmin, Rmax, Rmin };
        };

        // 4b. Factorize reduced stiffness matrix once into LU for O(n^2) solve per position
        const freeMap: number[] = [];
        for (let i = 0; i < nDOF; i++) {
            if (!constrained[i]) freeMap.push(i);
        }
        const nFree = freeMap.length;
        if (nFree < 1) throw new Error('The model has no free degrees of freedom.');
        const LU: number[][] = Array(nFree).fill(0).map(() => Array(nFree).fill(0));
        for (let i = 0; i < nFree; i++) {
            for (let j = 0; j < nFree; j++) {
                LU[i][j] = K_Global[freeMap[i]][freeMap[j]];
            }
        }
        // Doolittle LU without pivoting (reduced beam matrix is SPD). Guard zero pivots.
        for (let k = 0; k < nFree - 1; k++) {
            if (Math.abs(LU[k][k]) < 1e-12) throw new Error('Singular stiffness matrix - check span lengths and supports.');
            for (let i = k + 1; i < nFree; i++) {
                const factor = LU[i][k] / LU[k][k];
                LU[i][k] = factor;
                for (let j = k + 1; j < nFree; j++) {
                    LU[i][j] -= factor * LU[k][j];
                }
            }
        }
        if (Math.abs(LU[nFree - 1][nFree - 1]) < 1e-12) throw new Error('Singular stiffness matrix - check span lengths and supports.');

        const solveFactored = (F: number[]): number[] => {
            const y: number[] = Array(nFree);
            for (let i = 0; i < nFree; i++) y[i] = F[freeMap[i]];
            for (let i = 1; i < nFree; i++) {
                let sum = 0;
                for (let j = 0; j < i; j++) sum += LU[i][j] * y[j];
                y[i] -= sum;
            }
            for (let i = nFree - 1; i >= 0; i--) {
                let sum = 0;
                for (let j = i + 1; j < nFree; j++) sum += LU[i][j] * y[j];
                if (Math.abs(LU[i][i]) < 1e-12) throw new Error('Singular stiffness matrix during solve.');
                y[i] = (y[i] - sum) / LU[i][i];
            }
            const U_Full = Array(nDOF).fill(0);
            for (let i = 0; i < nFree; i++) U_Full[freeMap[i]] = y[i];
            return U_Full;
        };

        // 5-6. Single live-load case envelope (truck sweep + patterned UDL superposition).
        // VBA parity: Max uses truck+UDL_Max, Min uses truck+UDL_Min (not a single UDL value),
        // and reactions include UDL contributions.
        const runSingleCase = (factoredAxles: Axle[], w_udl: number) => {
            const udl = calculateUdlEnvelopes(w_udl);
            const shear = initEnvelope(nTotalElems * 2).map((p, i) => ({ ...p, x: xShear[i] }));
            const moment = initEnvelope(nNodes).map((p, i) => ({ ...p, x: xNodes[i] }));
            const deflect = initEnvelope(nNodes).map((p, i) => ({ ...p, x: xNodes[i] }));
            const react: ReactionEnvelope[] = supportPositions.map(x => ({ x, max: -Infinity, min: Infinity, govPos: 0 }));

            const fwdAxles = factoredAxles.map(a => ({ ...a }));
            // Deep-clone + reverse spacings (VBA revAxles parity). Never mutate fwdAxles via shared refs.
            const revAxles = fwdAxles.map(a => ({ ...a })).reverse();
            for (let i = 0; i < revAxles.length - 1; i++) {
                revAxles[i].spacing = fwdAxles[fwdAxles.length - 2 - i].spacing;
            }
            if (revAxles.length > 0) revAxles[revAxles.length - 1].spacing = 0;

            const truckLen = fwdAxles.reduce((acc, a) => acc + a.spacing, 0);
            const totalLen = xNodes[xNodes.length - 1];
            const startPos = -truckLen;
            const endPos = totalLen + truckLen;

            const runPass = (axles: Axle[]) => {
                const posSet = new Set<number>();
                for (let pos = startPos; pos <= endPos + 0.000001; pos += step) {
                    posSet.add(Math.round(pos * 100000) / 100000);
                }
                // Exact support alignments so axles land directly on every support
                for (const suppX of supportPositions) {
                    let dist = 0;
                    for (let k = 0; k < axles.length; k++) {
                        const lead = suppX + dist;
                        if (lead >= startPos - 0.001 && lead <= endPos + 0.001) {
                            posSet.add(Math.round(lead * 100000) / 100000);
                        }
                        if (k < axles.length - 1) dist += axles[k].spacing;
                    }
                }
                const sweepPositions = Array.from(posSet).sort((a, b) => a - b);

                for (const pos of sweepPositions) {
                    const F = Array(nDOF).fill(0);
                    const elemLoads = Array.from({ length: nTotalElems }, () => [0, 0, 0, 0]);

                    let axPos = pos;
                    this.applyPointLoad(F, elemLoads, axPos, axles[0].load, elemLens, xNodes);
                    for (let k = 0; k < axles.length - 1; k++) {
                        axPos -= axles[k].spacing;
                        this.applyPointLoad(F, elemLoads, axPos, axles[k + 1].load, elemLens, xNodes);
                    }

                    const U = solveFactored(F);
                    const { v, m, d, reactions } = this.calculateForces(U, nTotalElems, elemLens, E, I, K_Global, F, constrained, nElemsPerSpan, numSpans, elemLoads);

                    for (let i = 0; i < v.length; i++) {
                        const hi = v[i] + udl.Vmax[i];
                        if (hi > shear[i].max) shear[i].max = hi;
                        const lo = v[i] + udl.Vmin[i];
                        if (lo < shear[i].min) shear[i].min = lo;
                    }
                    for (let i = 0; i < m.length; i++) {
                        const hiM = m[i] + udl.Mmax[i];
                        if (hiM > moment[i].max) moment[i].max = hiM;
                        const loM = m[i] + udl.Mmin[i];
                        if (loM < moment[i].min) moment[i].min = loM;
                        const hiD = d[i] + udl.Dmax[i];
                        if (hiD > deflect[i].max) deflect[i].max = hiD;
                        const loD = d[i] + udl.Dmin[i];
                        if (loD < deflect[i].min) deflect[i].min = loD;
                    }
                    for (let i = 0; i < reactions.length; i++) {
                        const hiR = reactions[i] + udl.Rmax[i];
                        if (hiR > react[i].max) { react[i].max = hiR; react[i].govPos = pos; }
                        const loR = reactions[i] + udl.Rmin[i];
                        if (loR < react[i].min) react[i].min = loR;
                    }
                }
            };

            runPass(fwdAxles);
            runPass(revAxles);
            return { shear, moment, deflect, react };
        };

        // DLA applies to axle loads (including the 0.8 truck portion of lane load per
        // CSA S6-19 3.8.4.5 - truck portion carries DLA, UDL does not).
        const truckAxles = baseAxles.map(a => ({ ...a, load: a.load * (1 + dla) }));
        const laneAxles = baseAxles.map(a => ({ ...a, load: a.load * 0.8 * (1 + dla) }));
        const LANE_UDL = 9; // kN/m

        if (loadCase === 'truck') {
            const r = runSingleCase(truckAxles, 0);
            return { shear: r.shear, moment: r.moment, deflection: r.deflect, xNodes, reactions: r.react, supportPositions, dlaUsed: dla, incrementUsed: step, baseIncrement: truckIncrement };
        }
        if (loadCase === 'lane') {
            const r = runSingleCase(laneAxles, LANE_UDL);
            return { shear: r.shear, moment: r.moment, deflection: r.deflect, xNodes, reactions: r.react, supportPositions, dlaUsed: dla, incrementUsed: step, baseIncrement: truckIncrement };
        }
        // envelope: governing of truck-only and lane cases
        const t = runSingleCase(truckAxles, 0);
        const l = runSingleCase(laneAxles, LANE_UDL);
        const shear = t.shear.map((p, i) => ({ x: p.x, max: Math.max(p.max, l.shear[i].max), min: Math.min(p.min, l.shear[i].min) }));
        const moment = t.moment.map((p, i) => ({ x: p.x, max: Math.max(p.max, l.moment[i].max), min: Math.min(p.min, l.moment[i].min) }));
        const deflection = t.deflect.map((p, i) => ({ x: p.x, max: Math.max(p.max, l.deflect[i].max), min: Math.min(p.min, l.deflect[i].min) }));
        const reactions: ReactionEnvelope[] = t.react.map((p, i) => ({
            x: p.x,
            max: Math.max(p.max, l.react[i].max),
            min: Math.min(p.min, l.react[i].min),
            govPos: (l.react[i].max > p.max ? l.react[i].govPos : p.govPos) ?? 0,
        }));
        return { shear, moment, deflection, xNodes, reactions, supportPositions, dlaUsed: dla, incrementUsed: step, baseIncrement: truckIncrement };
    }

    private applyPointLoad(F: number[], elemLoads: number[][], pos: number, mag: number, elemLens: number[], xNodes: number[]) {
        const totalLen = xNodes[xNodes.length - 1];
        if (pos < -0.0001 || pos > totalLen + 0.0001) return;
        if (pos < 0) pos = 0;
        if (pos > totalLen) pos = totalLen;

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

        const f0 = -mag * N1 * 1000;
        const f1 = -mag * N2 * 1000;
        const f2 = -mag * N3 * 1000;
        const f3 = -mag * N4 * 1000;

        F[dofI] += f0;
        F[dofI + 1] += f1;
        F[dofJ] += f2;
        F[dofJ + 1] += f3;

        if (elemLoads && elemLoads[elemIdx]) {
            elemLoads[elemIdx][0] += f0;
            elemLoads[elemIdx][1] += f1;
            elemLoads[elemIdx][2] += f2;
            elemLoads[elemIdx][3] += f3;
        }
    }

    private calculateForces(
        U: number[],
        nElems: number,
        elemLens: number[],
        E: number,
        I: number,
        K_Global: number[][],
        F: number[],
        _constrained: boolean[],
        nElemsPerSpan: number,
        numSpans: number,
        elemLoads?: number[][]
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
                if (elemLoads && elemLoads[e]) {
                    f_loc[r] -= elemLoads[e][r];
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

        // Calculate reaction at each support node. R = K*U - F is already upward
        // positive (verified: central 100kN on 20m SSB gives +50kN each). VBA parity.
        for (const node of supportNodes) {
            const dof = node * 2; // Vertical DOF
            let reaction = 0;
            for (let j = 0; j < nDOF; j++) {
                reaction += K_Global[dof][j] * U[j];
            }
            reaction -= F[dof];
            reactions.push(reaction / 1000); // N -> kN, upward positive
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
    const [hoverData, setHoverData] = useState<{
        x: number;
        svgX: number;
        max: number;
        min: number;
    } | null>(null);

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

    // Find global extrema
    let globalMaxVal = -Infinity;
    let globalMaxX = xMin;
    let globalMinVal = Infinity;
    let globalMinX = xMin;

    for (let i = 0; i < data.length; i++) {
        const pt = data[i];
        const curMax = pt[dataKeyMax];
        const curMin = pt[dataKeyMin];
        if (curMax > globalMaxVal) {
            globalMaxVal = curMax;
            globalMaxX = pt.x;
        }
        if (curMin < globalMinVal) {
            globalMinVal = curMin;
            globalMinX = pt.x;
        }
    }

    const yRange = yMax - yMin;
    if (yRange === 0) { yMax += 1; yMin -= 1; }
    else { yMax += yRange * 0.1; yMin -= yRange * 0.1; }

    const xTicks = calculateTicks(xMin, xMax, 8);
    const yTicks = calculateTicks(yMin, yMax, 6);

    const xScale = (val: number) => padding.left + ((val - xMin) / (xMax - xMin)) * (width - padding.left - padding.right);
    const yScale = (val: number) =>
        flipY
            ? padding.top + ((val - yMin) / (yMax - yMin)) * (height - padding.top - padding.bottom)
            : height - padding.bottom - ((val - yMin) / (yMax - yMin)) * (height - padding.top - padding.bottom);

    const createPath = (vals: number[]) => {
        return vals.map((y, i) => `${i === 0 ? 'M' : 'L'} ${xScale(data[i].x)} ${yScale(y)}`).join(' ');
    };

    const pathMax = createPath(maxVals);
    const pathMin = createPath(minVals);

    const pathFill = `${pathMax} L ${xScale(data[data.length - 1].x)} ${yScale(minVals[minVals.length - 1])} ` +
        minVals.slice().reverse().map((y, i) => `L ${xScale(data[data.length - 1 - i].x)} ${yScale(y)}`).join(' ') + " Z";

    const zeroY = yScale(0);

    const unitSymbol = unit.includes('(') ? unit.split('(')[1].replace(')', '') : unit;

    const formatVal = (val: number) => {
        if (!Number.isFinite(val)) return '0.00';
        if (Math.abs(val) < 1e-6) return '0.00';
        if (Math.abs(val) < 0.005) return val.toFixed(4);
        return val.toFixed(2);
    };

    const formatValWithDetail = (val: number) => {
        if (!Number.isFinite(val)) return '0.00 ' + unitSymbol;
        if (unitSymbol.toLowerCase() === 'm') {
            const mm = (val * 1000).toFixed(2);
            return `${val.toFixed(4)} m (${mm} mm)`;
        }
        return `${formatVal(val)} ${unitSymbol}`;
    };

    const handlePointerMove = (e: React.PointerEvent<SVGSVGElement>) => {
        const svgRect = e.currentTarget.getBoundingClientRect();
        if (svgRect.width <= 0) return;
        const clientX = e.clientX - svgRect.left;
        const plotWidth = width - padding.left - padding.right;
        if (plotWidth <= 0) return;

        const rawRatio = (clientX - padding.left) / plotWidth;
        const clampedRatio = Math.max(0, Math.min(1, rawRatio));
        const targetX = xMin + clampedRatio * (xMax - xMin);

        // Find interpolated values or nearest segment
        let bestMax = 0;
        let bestMin = 0;
        let found = false;

        for (let i = 0; i < data.length - 1; i++) {
            const x1 = data[i].x;
            const x2 = data[i + 1].x;
            if (x1 === x2) continue; // vertical step boundary in stepped shear
            const minSegX = Math.min(x1, x2);
            const maxSegX = Math.max(x1, x2);
            if (targetX >= minSegX && targetX <= maxSegX) {
                const t = (targetX - x1) / (x2 - x1);
                bestMax = data[i][dataKeyMax] + t * (data[i + 1][dataKeyMax] - data[i][dataKeyMax]);
                bestMin = data[i][dataKeyMin] + t * (data[i + 1][dataKeyMin] - data[i][dataKeyMin]);
                found = true;
                break;
            }
        }

        if (!found) {
            let nearestDist = Infinity;
            let nearestIdx = 0;
            for (let i = 0; i < data.length; i++) {
                const dist = Math.abs(data[i].x - targetX);
                if (dist < nearestDist) {
                    nearestDist = dist;
                    nearestIdx = i;
                }
            }
            bestMax = data[nearestIdx][dataKeyMax];
            bestMin = data[nearestIdx][dataKeyMin];
        }

        const currentSvgX = xScale(targetX);
        setHoverData({
            x: targetX,
            svgX: currentSvgX,
            max: bestMax,
            min: bestMin
        });
    };

    const handlePointerLeave = () => {
        setHoverData(null);
    };

    const tooltipWidth = 205;
    const tooltipHeight = 68;
    const tooltipX = hoverData
        ? hoverData.svgX > width - tooltipWidth - 20
            ? Math.max(padding.left + 5, hoverData.svgX - tooltipWidth - 12)
            : Math.min(width - padding.right - tooltipWidth - 5, hoverData.svgX + 12)
        : 0;
    const tooltipY = Math.max(
        padding.top + 6,
        Math.min(height - padding.bottom - tooltipHeight - 6, padding.top + 8)
    );

    return (
        <div ref={containerRef} className="w-full bg-white rounded-lg shadow-sm border border-gray-200 p-4 mb-6" >
            <div className="flex flex-wrap items-center justify-between gap-2 mb-2">
                <h3 className="text-lg font-semibold text-gray-800">{title}</h3>
                <div className="flex flex-wrap items-center gap-2 text-xs">
                    <div className="flex items-center gap-1.5 px-2.5 py-1 rounded-md bg-blue-50 border border-blue-200 text-blue-800">
                        <span className="font-medium text-gray-500">Max:</span>
                        <span className="font-bold font-mono">{formatVal(globalMaxVal)} {unitSymbol}</span>
                        <span className="text-gray-500">@ {globalMaxX.toFixed(2)}m</span>
                    </div>
                    <div className="flex items-center gap-1.5 px-2.5 py-1 rounded-md bg-red-50 border border-red-200 text-red-800">
                        <span className="font-medium text-gray-500">Min:</span>
                        <span className="font-bold font-mono">{formatVal(globalMinVal)} {unitSymbol}</span>
                        <span className="text-gray-500">@ {globalMinX.toFixed(2)}m</span>
                    </div>
                </div>
            </div>

            <svg
                width={width}
                height={height}
                className="overflow-visible select-none"
                onPointerMove={handlePointerMove}
                onPointerLeave={handlePointerLeave}
            >

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

                <text
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

                <path d={pathFill} fill={color} fillOpacity="0.1" />
                <path d={pathMax} fill="none" stroke={color} strokeWidth="2" />
                <path d={pathMin} fill="none" stroke="#ef4444" strokeWidth="2" strokeDasharray="4 2" />

                {/* Extrema markers on the curves */}
                {Number.isFinite(globalMaxVal) && (
                    <circle
                        cx={xScale(globalMaxX)}
                        cy={yScale(globalMaxVal)}
                        r="4.5"
                        fill={color}
                        stroke="#ffffff"
                        strokeWidth="1.5"
                    >
                        <title>{`Max: ${formatVal(globalMaxVal)} ${unitSymbol} at x = ${globalMaxX.toFixed(2)} m`}</title>
                    </circle>
                )}
                {Number.isFinite(globalMinVal) && (
                    <circle
                        cx={xScale(globalMinX)}
                        cy={yScale(globalMinVal)}
                        r="4.5"
                        fill="#ef4444"
                        stroke="#ffffff"
                        strokeWidth="1.5"
                    >
                        <title>{`Min: ${formatVal(globalMinVal)} ${unitSymbol} at x = ${globalMinX.toFixed(2)} m`}</title>
                    </circle>
                )}

                {/* Dynamic Cursor Overlay */}
                {hoverData && (
                    <g pointerEvents="none">
                        <line
                            x1={hoverData.svgX}
                            y1={padding.top}
                            x2={hoverData.svgX}
                            y2={height - padding.bottom}
                            stroke="#475569"
                            strokeWidth="1.2"
                            strokeDasharray="3 3"
                        />
                        <circle
                            cx={hoverData.svgX}
                            cy={yScale(hoverData.max)}
                            r="4.5"
                            fill={color}
                            stroke="#ffffff"
                            strokeWidth="2"
                        />
                        <circle
                            cx={hoverData.svgX}
                            cy={yScale(hoverData.min)}
                            r="4.5"
                            fill="#ef4444"
                            stroke="#ffffff"
                            strokeWidth="2"
                        />

                        {/* Tooltip Card */}
                        <rect
                            x={tooltipX}
                            y={tooltipY}
                            width={tooltipWidth}
                            height={tooltipHeight}
                            rx="6"
                            fill="#0f172a"
                            opacity="0.94"
                        />
                        <text x={tooltipX + 10} y={tooltipY + 18} fill="#f8fafc" fontSize="11" fontWeight="700">
                            x = {hoverData.x.toFixed(2)} m
                        </text>
                        <text x={tooltipX + 10} y={tooltipY + 36} fill="#60a5fa" fontSize="10" fontWeight="500">
                            Max: {formatValWithDetail(hoverData.max)}
                        </text>
                        <text x={tooltipX + 10} y={tooltipY + 54} fill="#f87171" fontSize="10" fontWeight="500">
                            Min: {formatValWithDetail(hoverData.min)}
                        </text>
                    </g>
                )}

                {/* Interactive Overlay Layer */}
                <rect
                    x={padding.left}
                    y={padding.top}
                    width={width - padding.left - padding.right}
                    height={height - padding.top - padding.bottom}
                    fill="transparent"
                    className="cursor-crosshair"
                />

            </svg>
            <div className="flex justify-center gap-6 mt-2 text-sm" >
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
    supportPositions,
    dla
}: {
    spans: Span[],
    reactions: ReactionEnvelope[],
    supportPositions: number[],
    dla?: number
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
        // Reactions are upward-positive kN (sign fixed). Show governing max; min (uplift) in title.
        const displayR = reaction.max;

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
                    <title>{`Max ${reaction.max.toFixed(1)} kN, Min ${reaction.min.toFixed(1)} kN${reaction.govPos !== undefined ? `, gov at ${reaction.govPos.toFixed(2)}m` : ''}`}</title>
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
                <span>↑ Maximum reaction forces shown{dla !== undefined ? ` (Applied DLA = ${(dla * 100).toFixed(0)}%)` : ''}</span>
            </div>
        </div>
    );
};

// Helper to calculate automated DLA based on CSA S6-19 Cl. 3.8.4.5 and span arrangement.
// Conservative: evaluates each span and governs with the highest DLA (short spans control).
export function computeAutoDlaInfo(spans: Span[], axles: Axle[]): {
    dla: number;
    axleCount: number;
    maxSpan: number;
    desc: string;
} {
    const maxSpan = Math.max(...spans.map(s => s.length), 0);
    const countOn = (len: number) => {
        let best = 1;
        for (let i = 0; i < axles.length; i++) {
            let count = 1;
            let cum = 0;
            for (let j = i; j < axles.length - 1; j++) {
                cum += axles[j].spacing;
                if (cum <= len + 0.000001) count++;
                else break;
            }
            if (count > best) best = count;
        }
        return best;
    };
    const dlaFor = (n: number) => (n <= 1 ? 0.40 : n === 2 ? 0.30 : 0.25);
    let dla = 0.25;
    let govCount = 3;
    for (const s of spans) {
        const n = countOn(s.length);
        const d = dlaFor(n);
        if (d > dla) {
            dla = d;
            govCount = n;
        }
    }
    let desc = '≥ 3 axles on span';
    if (govCount <= 1) desc = '1 axle on span';
    else if (govCount === 2) desc = '2 axles on span (tandem)';
    return { dla, axleCount: govCount, maxSpan, desc };
}

export default function BeamAnalysisApp() {
    const [spans, setSpans] = useState<Span[]>(DEFAULT_SPANS);
    const [axles, _setAxles] = useState<Axle[]>(DEFAULT_AXLES);
    const [config, setConfig] = useState<AnalysisConfig>(DEFAULT_CONFIG);
    const [results, setResults] = useState<AnalysisResults | null>(null);
    const [isAnalyzing, setIsAnalyzing] = useState(false);
    const [activeTab, setActiveTab] = useState<'config' | 'results'>('config');

    const autoDlaInfo = useMemo(() => computeAutoDlaInfo(spans, axles), [spans, axles]);
    const stepInfo = useMemo(
        () => computeEffectiveIncrement(spans, axles, config.truckIncrement, config.nElemsPerSpan),
        [spans, axles, config.truckIncrement, config.nElemsPerSpan]
    );

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
                alert(`Analysis failed: ${e instanceof Error ? e.message : 'Check inputs.'}`);
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

        // Reactions Sheet (VBA parity: summary with gov truck position)
        const reactionData = results.reactions.map((r, i) => ({
            "Support": `Support ${i + 1}`,
            "Location (m)": parseFloat(r.x.toFixed(3)),
            "Max Reaction (kN)": parseFloat(r.max.toFixed(3)),
            "Min Reaction (kN)": parseFloat(r.min.toFixed(3)),
            "Gov Truck Pos (m)": r.govPos !== undefined ? parseFloat(r.govPos.toFixed(3)) : 0,
        }));
        // @ts-ignore
        const wsReact = window.XLSX.utils.json_to_sheet(reactionData);
        // @ts-ignore
        window.XLSX.utils.book_append_sheet(wb, wsReact, "Reactions");

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
                                                    onChange={(e) => setConfig({ ...config, loadCase: e.target.value as 'truck' | 'lane' | 'envelope' })}
                                                    className="w-full bg-white border border-gray-300 rounded px-3 py-2 text-sm focus:ring-2 focus:ring-blue-500 outline-none"
                                                >
                                                    <option value="truck" > CL-625 Truck Only (Standard) </option>
                                                    <option value="lane" > CL-625 Lane Load (80% Truck with DLA + 9 kN/m patterned) </option>
                                                    <option value="envelope" > Envelope (max of Truck and Lane) </option>
                                                </select>
                                            </div>

                                            <div>
                                                <div className="flex items-center justify-between mb-1.5">
                                                    <label className="block text-sm font-medium text-gray-700">
                                                        Dynamic Load Allowance (DLA)
                                                    </label>
                                                    <label className="flex items-center gap-1.5 text-xs text-gray-600 cursor-pointer select-none">
                                                        <input
                                                            type="checkbox"
                                                            checked={config.dlaOverride !== null && config.dlaOverride !== undefined}
                                                            onChange={(e) => {
                                                                const checked = e.target.checked;
                                                                setConfig({
                                                                    ...config,
                                                                    dlaOverride: checked ? (config.dlaOverride ?? autoDlaInfo.dla) : null,
                                                                });
                                                            }}
                                                            className="rounded border-gray-300 text-blue-600 focus:ring-blue-500 h-3.5 w-3.5 cursor-pointer"
                                                        />
                                                        <span className="font-medium">Override DLA</span>
                                                    </label>
                                                </div>

                                                {config.dlaOverride === null || config.dlaOverride === undefined ? (
                                                    <div className="bg-slate-50 border border-slate-200 rounded p-2.5 flex items-center justify-between">
                                                        <div className="flex items-center gap-2">
                                                            <span className="inline-flex items-center px-2 py-0.5 rounded text-xs font-bold bg-blue-100 text-blue-800">
                                                                Auto: {(autoDlaInfo.dla * 100).toFixed(0)}% (DLA = {autoDlaInfo.dla.toFixed(2)})
                                                            </span>
                                                            <span className="text-xs text-slate-600">
                                                                {autoDlaInfo.desc} (L<sub>max</sub> = {autoDlaInfo.maxSpan.toFixed(2)}m)
                                                            </span>
                                                        </div>
                                                        <span className="text-[10px] text-slate-400 font-medium">CSA S6 Cl. 3.8.4.5</span>
                                                    </div>
                                                ) : (
                                                    <div className="space-y-1.5">
                                                        <div className="flex items-center gap-2">
                                                            <input
                                                                type="number"
                                                                step="0.01"
                                                                min="0"
                                                                max="1.0"
                                                                value={config.dlaOverride}
                                                                onChange={(e) => {
                                                                    const val = parseFloat(e.target.value);
                                                                    setConfig({
                                                                        ...config,
                                                                        dlaOverride: isNaN(val) ? 0 : val,
                                                                    });
                                                                }}
                                                                className="w-1/3 bg-white border border-blue-400 rounded px-2.5 py-1 text-sm focus:ring-2 focus:ring-blue-500 outline-none font-semibold text-blue-900"
                                                                placeholder="e.g. 0.30"
                                                            />
                                                            <span className="text-xs text-amber-800 bg-amber-50 border border-amber-200 rounded px-2.5 py-1">
                                                                Manual override active: <strong>{((config.dlaOverride ?? 0) * 100).toFixed(1)}%</strong>
                                                            </span>
                                                        </div>
                                                    </div>
                                                )}
                                                <span className="text-[11px] text-gray-500 mt-1 block">
                                                    Automated per CSA S6 Cl. 3.8.4.5 based on span arrangement and axles on span (40% for 1 axle, 30% for 2-axle tandem, 25% for ≥ 3 axles).
                                                </span>
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

                                            <div className="grid grid-cols-2 gap-4" >
                                                <div>
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Elements per span </label>
                                                    < input
                                                        type="number"
                                                        min="2"
                                                        max="200"
                                                        step="1"
                                                        value={config.nElemsPerSpan}
                                                        onChange={(e) => setConfig({ ...config, nElemsPerSpan: Math.max(2, Math.round(parseFloat(e.target.value) || 40)) })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                                <div>
                                                    <label className="block text-sm font-medium text-gray-700 mb-1" > Truck step, base (m) </label>
                                                    < input
                                                        type="number"
                                                        min="0.02"
                                                        max="2"
                                                        step="0.05"
                                                        value={config.truckIncrement}
                                                        onChange={(e) => setConfig({ ...config, truckIncrement: parseFloat(e.target.value) || 0.25 })}
                                                        className="w-full border border-gray-300 rounded px-2 py-1 text-sm"
                                                    />
                                                </div>
                                            </div>

                                            < div className="bg-blue-50 p-3 rounded text-sm text-blue-800 flex gap-2 items-start" >
                                                <AlertCircle className="w-4 h-4 mt-0.5 shrink-0" />
                                                <p>Mesh: {config.nElemsPerSpan} elements per span. Truck sweep uses <strong>{stepInfo.effective.toFixed(3)}m</strong> steps{stepInfo.wasReduced ? <> (auto-refined from {config.truckIncrement}m — {stepInfo.reason})</> : <> (base {config.truckIncrement}m)</>}. Support alignments are always added exactly.</p>
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
                            {results.incrementUsed !== undefined && (
                                <div className="w-full bg-blue-50 border border-blue-200 rounded-lg px-4 py-2 mb-4 text-xs text-blue-800">
                                    Sweep step used: <strong>{results.incrementUsed.toFixed(3)}m</strong>
                                    {results.baseIncrement !== undefined && Math.abs(results.incrementUsed - results.baseIncrement) > 1e-9
                                        ? <> (auto-refined from base {results.baseIncrement.toFixed(3)}m for short spans/axle spacing)</>
                                        : <> (base setting)</>}.
                                    Support peaks captured exactly via alignment positions.
                                </div>
                            )}
                            <BeamReactionDiagram
                                spans={spans}
                                reactions={results.reactions}
                                supportPositions={results.supportPositions}
                                dla={results.dlaUsed}
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
                                    Download the analysis results as an Excel (.xlsx) file with separate sheets for Shear, Moment, Deflection, and Support Reactions.
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