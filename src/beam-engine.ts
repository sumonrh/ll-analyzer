export type Span = { id: string; length: number };
export type Axle = { id: string; load: number; spacing: number };
export type LoadCase = 'truck' | 'lane' | 'envelope';
export type AnalysisConfig = {
    E: number;
    I: number;
    nElemsPerSpan: number;
    truckIncrement: number;
    loadCase: LoadCase;
    dlaOverride?: number | null;
    dlaMultiplier?: number;
    laneUdl?: number | null;
};
export type EnvelopePoint = { x: number; max: number; min: number };
export type ReactionEnvelope = EnvelopePoint & { govPos: number };
export type UdlInterval = {
    envelope: 'max' | 'min';
    element: number;
    xiStart: number;
    xiEnd: number;
    start: number;
    end: number;
    contribution: number;
};
export type UdlTracer = {
    response: 'moment';
    nodeIndex: number;
    x: number;
    unit: 'kNm';
    intensity: number;
    max: number;
    min: number;
    reconstructedMax: number;
    reconstructedMin: number;
    tolerance: number;
    intervals: UdlInterval[];
    influenceLine: (EnvelopePoint & { ordinate: number })[];
};
export type CaseResults = {
    shear: EnvelopePoint[];
    moment: EnvelopePoint[];
    deflection: EnvelopePoint[];
    reactions: ReactionEnvelope[];
    reactionDiagrams: EnvelopePoint[][];
    dlaUsed: number;
    dlaAuto: boolean;
    dlaBase: number;
    dlaMultiplier: number;
};
export type AnalysisResults = CaseResults & {
    cases: Partial<Record<LoadCase, CaseResults>>;
    loadCase: LoadCase;
    spans: Span[];
    axles: Axle[];
    config: AnalysisConfig;
    xNodes: number[];
    supportPositions: number[];
    truckPositions: number[];
    incrementUsed: number;
    baseIncrement: number;
    incrementReason: string;
    elapsedMs: number;
    udlTracer?: UdlTracer;
    stats: { factorizations: number; truckSolves: number; influenceSolves: number; udlIntegrations: number };
};
export type AnalysisProgress = { fraction: number; message: string };
export type AnalysisRequest = { spans: Span[]; axles: Axle[]; config: AnalysisConfig };
export type AnalysisResponse =
    | { type: 'progress'; progress: AnalysisProgress }
    | { type: 'result'; result: AnalysisResults }
    | { type: 'error'; message: string };

export const DEFAULT_SPANS: Span[] = [
    { id: 's1', length: 20 }, { id: 's2', length: 25 }, { id: 's3', length: 20 },
];
export const DEFAULT_AXLES: Axle[] = [
    { id: 'a1', load: 50, spacing: 3.6 },
    { id: 'a2', load: 125, spacing: 1.2 },
    { id: 'a3', load: 125, spacing: 6.6 },
    { id: 'a4', load: 175, spacing: 6.6 },
    { id: 'a5', load: 150, spacing: 0 },
];
export const DEFAULT_CONFIG: AnalysisConfig = {
    E: 200000000000, I: 0.005, nElemsPerSpan: 40, truckIncrement: 0.25,
    loadCase: 'truck', dlaOverride: null, dlaMultiplier: 1, laneUdl: 9,
};
export const MAX_PATTERN_SPANS = 12;
export const MAX_AXLES = 20;
export const MAX_SWEEP_STEPS = 6000;
const MIN_SWEEP_STEP = 0.02;
const BAND = 3;
// Load rules from the supplied VBA reference; MIDAS parity needs a matching benchmark.
const DLA_FACTOR = 0.25;
const LANE_TRUCK_FACTOR = 0.8;
const LANE_UDL = 9;
const MAX_TOTAL_ELEMENTS = 1000;
const INFLUENCE_ROOT_TOL = 1e-10;
const INFLUENCE_VALUE_TOL = 1e-12;
const UDL_VERIFY_TOL = 1e-7;

function maxAxlesOnLength(length: number, axles: Axle[]): number {
    let best = 1;
    for (let first = 0; first < axles.length; first++) {
        let distance = 0;
        let count = 1;
        for (let next = first; next < axles.length - 1; next++) {
            distance += axles[next].spacing;
            if (distance > length + 1e-6) break;
            count++;
        }
        best = Math.max(best, count);
    }
    return best;
}

export function computeAutoDlaInfo(spans: Span[], axles: Axle[]) {
    let dla = 0.25;
    let axleCount = 3;
    let governingSpan = spans[0]?.length ?? 0;
    for (const span of spans) {
        const count = maxAxlesOnLength(span.length, axles);
        const candidate = count <= 1 ? 0.40 : count === 2 ? 0.30 : 0.25;
        if (candidate > dla) {
            dla = candidate;
            axleCount = count;
            governingSpan = span.length;
        }
    }
    const desc = axleCount <= 1 ? '1 axle on span'
        : axleCount === 2 ? '2 axles on span (tandem)' : '>= 3 axles on span';
    return { dla, axleCount, maxSpan: Math.max(0, ...spans.map(s => s.length)), governingSpan, desc };
}

export function resolveDla(config: AnalysisConfig): { isAuto: boolean; base: number; multiplier: number; effective: number } {
    let multiplier = config.dlaMultiplier ?? 1;
    if (!Number.isFinite(multiplier)) multiplier = 1;
    multiplier = Math.max(0, Math.min(1, multiplier));
    const isAuto = config.dlaOverride == null;
    const base = isAuto ? 0 : (config.dlaOverride as number);
    return { isAuto, base, multiplier, effective: base * multiplier };
}

export function truckGroupDla(count: number, frontThree: boolean): number {
    if (count <= 0) return 0;
    if (count === 1) return 0.40;
    if (count === 2) return 0.30;
    if (count === 3) return frontThree ? 0.30 : DLA_FACTOR;
    return DLA_FACTOR;
}

export function computeEffectiveIncrement(
    spans: Span[], axles: Axle[], baseIncrement: number, nElemsPerSpan: number
) {
    if (!Number.isFinite(baseIncrement) || baseIncrement <= 0 ||
        !Number.isInteger(nElemsPerSpan) || nElemsPerSpan < 2 ||
        spans.length === 0 || spans.some(s => !Number.isFinite(s.length) || s.length <= 0)) {
        return { effective: baseIncrement, reason: 'Enter valid analysis settings.', wasReduced: false, wasAdjusted: false };
    }
    const minSpan = Math.min(...spans.map(s => s.length));
    const gaps = axles.slice(0, -1).map(a => a.spacing).filter(g => g > 0 && Number.isFinite(g));
    const truckLength = axles.slice(0, -1).reduce((sum, a) => sum + a.spacing, 0);
    const sweepLength = spans.reduce((sum, s) => sum + s.length, 0) + 2 * truckLength;
    let effective = baseIncrement;
    let reason = 'base setting';
    for (const [candidate, label] of [
        [minSpan / 40, 'short span control'],
        [minSpan / nElemsPerSpan, 'element control'],
        [gaps.length ? Math.min(...gaps) / 8 : Infinity, 'axle-spacing control'],
    ] as const) {
        if (candidate < effective) {
            effective = candidate;
            reason = label;
        }
    }
    if (effective < MIN_SWEEP_STEP) {
        effective = MIN_SWEEP_STEP;
        reason = 'minimum step floor (0.02m)';
    }
    if (sweepLength / effective > MAX_SWEEP_STEPS) {
        effective = sweepLength / MAX_SWEEP_STEPS;
        reason = `step-count cap (${MAX_SWEEP_STEPS} intervals per direction)`;
    }
    const safeMinimum = Math.max(MIN_SWEEP_STEP, sweepLength / MAX_SWEEP_STEPS);
    effective = Math.ceil(effective * 200 - 1e-10) / 200;
    if (baseIncrement >= safeMinimum) effective = Math.min(effective, baseIncrement);
    return {
        effective, reason,
        wasReduced: effective < baseIncrement - 1e-9,
        wasAdjusted: Math.abs(effective - baseIncrement) > 1e-9,
    };
}

export function validateInputs({ spans, axles, config }: AnalysisRequest): void {
    if (spans.length < 1) throw new Error('At least one span is required.');
    spans.forEach((span, i) => {
        if (!Number.isFinite(span.length) || span.length <= 0)
            throw new Error(`Span ${i + 1} length must be greater than zero.`);
    });
    if (!Number.isFinite(config.E) || config.E <= 0)
        throw new Error('Elastic modulus E must be greater than zero.');
    if (!Number.isFinite(config.I) || config.I <= 0)
        throw new Error('Moment of inertia I must be greater than zero.');
    if (!Number.isInteger(config.nElemsPerSpan) || config.nElemsPerSpan < 2)
        throw new Error('Elements per span must be an integer >= 2.');
    if (!Number.isFinite(config.truckIncrement) || config.truckIncrement <= 0)
        throw new Error('Truck increment must be greater than zero.');
    if (!['truck', 'lane', 'envelope'].includes(config.loadCase))
        throw new Error('Choose Truck, Lane Load or Envelope.');
    if (config.loadCase !== 'truck' && spans.length > MAX_PATTERN_SPANS)
        throw new Error(`Lane load patterning is limited to ${MAX_PATTERN_SPANS} spans.`);
    if (axles.length < 1 || axles.length > MAX_AXLES)
        throw new Error(`Enter between 1 and ${MAX_AXLES} axles.`);
    axles.forEach((axle, i) => {
        if (!Number.isFinite(axle.load) || axle.load < 0)
            throw new Error(`Axle ${i + 1} load must be >= 0.`);
        if (i < axles.length - 1 && (!Number.isFinite(axle.spacing) || axle.spacing < 0))
            throw new Error(`Axle ${i + 1} spacing must be >= 0.`);
    });
    if (config.dlaOverride != null && (!Number.isFinite(config.dlaOverride) || config.dlaOverride < 0))
        throw new Error('DLA must be >= 0.');
    const mult = config.dlaMultiplier ?? 1;
    if (!Number.isFinite(mult) || mult < 0 || mult > 1)
        throw new Error('Truck DLA multiplier (d) must be between 0 and 1.');
    const laneUdl = config.laneUdl ?? LANE_UDL;
    if (!Number.isFinite(laneUdl) || laneUdl < 0)
        throw new Error('Lane UDL must be >= 0 kN/m.');
    const totalLength = spans.reduce((sum, s) => sum + s.length, 0);
    const truckLength = axles.slice(0, -1).reduce((sum, a) => sum + a.spacing, 0);
    if (!Number.isFinite(totalLength + 2 * truckLength) || !Number.isFinite(config.E * config.I))
        throw new Error('Geometry or stiffness exceeds the supported numeric range.');
    if (!Number.isSafeInteger(spans.length * config.nElemsPerSpan) ||
        spans.length * config.nElemsPerSpan > MAX_TOTAL_ELEMENTS)
        throw new Error(`The model is limited to ${MAX_TOTAL_ELEMENTS} beam elements. Reduce the mesh resolution.`);
}

function reverseAxles(axles: Axle[]): Axle[] {
    return axles.map((_, i) => ({
        ...axles[axles.length - 1 - i],
        spacing: i < axles.length - 1 ? axles[axles.length - 2 - i].spacing : 0,
    }));
}

export function buildTruckPositions(supports: number[], axles: Axle[], step: number): number[] {
    const truckLength = axles.slice(0, -1).reduce((sum, a) => sum + a.spacing, 0);
    const start = -truckLength;
    const end = supports[supports.length - 1] + truckLength;
    const count = Math.max(2, Math.floor((end - start) / step + 1e-6) + 1);
    const uniform = Array.from({ length: count }, (_, i) => Math.min(end, start + i * step));
    uniform[uniform.length - 1] = end;
    const positions = new Set(uniform);
    for (const direction of [axles, reverseAxles(axles)]) {
        for (const support of supports) {
            let distance = 0;
            for (let i = 0; i < direction.length; i++) {
                for (const sign of [1, -1]) {
                    const lead = support + sign * distance;
                    if (lead >= start && lead <= end) positions.add(lead);
                }
                if (i < direction.length - 1) distance += direction[i].spacing;
            }
        }
    }
    // Keep alignment coordinates exact; do not round them into a nearby reaction bin.
    return [...positions].sort((a, b) => a - b)
        .filter((value, i, all) => i === 0 || value - all[i - 1] > 1e-9);
}

class BeamSystem {
    readonly xNodes: number[] = [0];
    readonly xShear: number[] = [];
    readonly supports: number[] = [0];
    readonly nElems: number;
    readonly nNodes: number;
    readonly nDOF: number;
    readonly momentOffset: number;
    readonly deflectionOffset: number;
    readonly reactionOffset: number;
    readonly response: Float64Array;
    readonly lengths: Float64Array;
    private readonly localK: Float64Array;
    private readonly stiffness: Float64Array;
    private readonly freeMap: number[] = [];
    private readonly lower: Float64Array;
    private readonly diagonal: Float64Array;
    private readonly load: Float64Array;
    private readonly elemLoads: Float64Array;
    private readonly displacement: Float64Array;
    private readonly work: Float64Array;
    private readonly endForces: Float64Array;
    private readonly supportDOFs: number[] = [0];
    private readonly config: AnalysisConfig;

    constructor(spans: Span[], config: AnalysisConfig) {
        this.config = config;
        this.nElems = spans.length * config.nElemsPerSpan;
        this.nNodes = this.nElems + 1;
        this.nDOF = this.nNodes * 2;
        this.momentOffset = this.nElems * 2;
        this.deflectionOffset = this.momentOffset + this.nNodes;
        this.reactionOffset = this.deflectionOffset + this.nNodes;
        this.response = new Float64Array(this.reactionOffset + spans.length + 1);
        this.lengths = new Float64Array(this.nElems);
        this.localK = new Float64Array(this.nElems * 16);
        this.stiffness = new Float64Array(this.nDOF * 4);
        this.load = new Float64Array(this.nDOF);
        this.elemLoads = new Float64Array(this.nElems * 4);
        this.displacement = new Float64Array(this.nDOF);
        this.endForces = new Float64Array(this.nElems * 4);
        let start = 0;
        let elem = 0;
        for (const span of spans) {
            const length = span.length / config.nElemsPerSpan;
            for (let i = 0; i < config.nElemsPerSpan; i++, elem++) {
                const x = start + (i + 1) * length;
                if (x <= this.xNodes[this.xNodes.length - 1])
                    throw new Error('Span lengths cannot be resolved at this mesh scale.');
                this.xShear.push(this.xNodes[this.xNodes.length - 1], x);
                this.xNodes.push(x);
                this.lengths[elem] = length;
                const l2 = length * length;
                const coeff = config.E * config.I / (length ** 3);
                const matrix = [
                    12, 6 * length, -12, 6 * length,
                    6 * length, 4 * l2, -6 * length, 2 * l2,
                    -12, -6 * length, 12, -6 * length,
                    6 * length, 2 * l2, -6 * length, 4 * l2,
                ];
                for (let r = 0; r < 4; r++) {
                    for (let c = 0; c < 4; c++) {
                        const value = matrix[r * 4 + c] * coeff;
                        if (!Number.isFinite(value))
                            throw new Error('Element stiffness exceeds the supported numeric range.');
                        this.localK[elem * 16 + r * 4 + c] = value;
                        if (r >= c) this.stiffness[(elem * 2 + r) * 4 + r - c] += value;
                    }
                }
            }
            start += span.length;
            this.supports.push(start);
            this.supportDOFs.push(elem * 2);
        }
        const constrained = new Set(this.supportDOFs);
        for (let i = 0; i < this.nDOF; i++) if (!constrained.has(i)) this.freeMap.push(i);
        this.lower = new Float64Array(this.freeMap.length * 4);
        this.diagonal = new Float64Array(this.freeMap.length);
        this.work = new Float64Array(this.freeMap.length);
        this.factorize();
    }

    private k(row: number, col: number): number {
        const hi = Math.max(row, col);
        const gap = Math.abs(row - col);
        return gap <= BAND ? this.stiffness[hi * 4 + gap] : 0;
    }

    private factorize(): void {
        // The reduced Euler-Bernoulli matrix is SPD with half-bandwidth <= 3.
        // Banded LDL^T is equivalent to the VBA cached LU, without dense zero work.
        for (let i = 0; i < this.freeMap.length; i++) {
            for (let j = Math.max(0, i - BAND); j < i; j++) {
                let value = this.k(this.freeMap[i], this.freeMap[j]);
                for (let k = Math.max(0, i - BAND); k < j; k++)
                    value -= this.lower[i * 4 + i - k] * this.diagonal[k] * this.lower[j * 4 + j - k];
                this.lower[i * 4 + i - j] = value / this.diagonal[j];
            }
            let pivot = this.k(this.freeMap[i], this.freeMap[i]);
            for (let j = Math.max(0, i - BAND); j < i; j++)
                pivot -= this.lower[i * 4 + i - j] ** 2 * this.diagonal[j];
            if (!Number.isFinite(pivot) || pivot <= 0)
                throw new Error('Singular or ill-conditioned stiffness matrix. Check span lengths and supports.');
            this.diagonal[i] = pivot;
        }
    }

    private clearLoads(): void {
        this.load.fill(0);
        this.elemLoads.fill(0);
    }

    private pointLoad(position: number, magnitude: number): void {
        const totalLength = this.supports[this.supports.length - 1];
        if (position < -1e-4 || position > totalLength + 1e-4) return;
        const pos = Math.max(0, Math.min(totalLength, position));
        let low = 0;
        let high = this.nElems - 1;
        while (low < high) {
            const mid = Math.floor((low + high) / 2);
            if (pos <= this.xNodes[mid + 1] + 1e-6) high = mid;
            else low = mid + 1;
        }
        const le = this.lengths[low];
        const xi = Math.max(0, Math.min(1, (pos - this.xNodes[low]) / le));
        this.elementPointLoad(low, xi, magnitude);
    }

    private elementPointLoad(element: number, xi: number, magnitude: number): void {
        const le = this.lengths[element];
        const x2 = xi * xi;
        const x3 = x2 * xi;
        const shapes = [1 - 3 * x2 + 2 * x3, le * (xi - 2 * x2 + x3),
            3 * x2 - 2 * x3, le * (-x2 + x3)];
        for (let i = 0; i < 4; i++) {
            const value = -magnitude * shapes[i] * 1000;
            this.load[element * 2 + i] += value;
            this.elemLoads[element * 4 + i] += value;
        }
    }

    truck(position: number, axles: Axle[]): Float64Array {
        this.clearLoads();
        let axlePosition = position;
        for (let i = 0; i < axles.length; i++) {
            this.pointLoad(axlePosition, axles[i].load);
            if (i < axles.length - 1) axlePosition -= axles[i].spacing;
        }
        return this.solve();
    }

    /** Solve a unit (1 kN) point load and return a copy of the packed response. */
    unitResponse(element: number, xi: number): Float64Array {
        this.clearLoads();
        this.elementPointLoad(element, xi, 1);
        this.solve();
        return this.response.slice();
    }

    partialUdlResponse(intervals: UdlInterval[], envelope: 'max' | 'min', w: number): Float64Array {
        this.clearLoads();
        for (const interval of intervals) {
            if (interval.envelope !== envelope) continue;
            const e = interval.element - 1;
            const le = this.lengths[e];
            const a = interval.xiStart;
            const b = interval.xiEnd;
            const i0 = b - a;
            const i1 = i0 * (a + b) / 2;
            const i2 = i0 * (a * a + a * b + b * b) / 3;
            const i3 = i0 * (a ** 3 + a * a * b + a * b * b + b ** 3) / 4;
            const factor = -w * le * 1000;
            const values = [
                factor * (i0 - 3 * i2 + 2 * i3),
                factor * le * (i1 - 2 * i2 + i3),
                factor * (3 * i2 - 2 * i3),
                factor * le * (-i2 + i3),
            ];
            for (let i = 0; i < 4; i++) {
                this.load[e * 2 + i] += values[i];
                this.elemLoads[e * 4 + i] += values[i];
            }
        }
        return this.solve();
    }

    private solve(): Float64Array {
        for (let i = 0; i < this.freeMap.length; i++) {
            let value = this.load[this.freeMap[i]];
            for (let j = Math.max(0, i - BAND); j < i; j++)
                value -= this.lower[i * 4 + i - j] * this.work[j];
            this.work[i] = value;
        }
        for (let i = this.freeMap.length - 1; i >= 0; i--) {
            let value = this.work[i] / this.diagonal[i];
            for (let j = i + 1; j <= Math.min(this.freeMap.length - 1, i + BAND); j++)
                value -= this.lower[j * 4 + j - i] * this.displacement[this.freeMap[j]];
            this.displacement[this.freeMap[i]] = value;
        }
        for (let e = 0; e < this.nElems; e++) {
            for (let r = 0; r < 4; r++) {
                let force = -this.elemLoads[e * 4 + r];
                for (let c = 0; c < 4; c++)
                    force += this.localK[e * 16 + r * 4 + c] * this.displacement[e * 2 + c];
                this.endForces[e * 4 + r] = force;
            }
            this.response[e * 2] = this.endForces[e * 4] / 1000;
            this.response[e * 2 + 1] = -this.endForces[e * 4 + 2] / 1000;
        }
        for (let n = 0; n < this.nNodes; n++) {
            const moment = n === 0 ? this.endForces[1]
                : n === this.nElems ? -this.endForces[(n - 1) * 4 + 3]
                    : (this.endForces[n * 4 + 1] - this.endForces[(n - 1) * 4 + 3]) / 2;
            this.response[this.momentOffset + n] = -moment / 1000;
            this.response[this.deflectionOffset + n] = this.displacement[n * 2];
        }
        this.supportDOFs.forEach((dof, s) => {
            let reaction = -this.load[dof];
            for (let j = Math.max(0, dof - BAND); j <= Math.min(this.nDOF - 1, dof + BAND); j++)
                reaction += this.k(dof, j) * this.displacement[j];
            this.response[this.reactionOffset + s] = reaction / 1000;
        });
        if (this.response.some(value => !Number.isFinite(value)))
            throw new Error('The analysis produced non-finite results. Check geometry, stiffness and loads.');
        return this.response;
    }
}

// ---------------------------------------------------------------------------
// VBA-parity influence engine
// FEA + influence-line UDL placement, continuous truck optimisation,
// placement-dependent selected-axle DLA with d multiplier.
// ---------------------------------------------------------------------------

function fitInfluenceCubic(y0: number, y1: number, y2: number, y3: number): number[] {
    const d1 = y1 - y0;
    const d2 = y2 - 2 * y1 + y0;
    const d3 = y3 - 3 * y2 + 3 * y1 - y0;
    return [
        y0 - 0.5 * d1 + 0.375 * d2 - 0.3125 * d3,
        4 * d1 - 4 * d2 + 23 * d3 / 6,
        8 * d2 - 12 * d3,
        32 * d3 / 3,
    ];
}

function influenceValue(coeff: ArrayLike<number>, xi: number): number {
    return ((coeff[3] * xi + coeff[2]) * xi + coeff[1]) * xi + coeff[0];
}

function influenceIntegral(coeff: ArrayLike<number>, a: number, b: number): number {
    return (b - a) * (coeff[0] + coeff[1] * (a + b) / 2 +
        coeff[2] * (a * a + a * b + b * b) / 3 +
        coeff[3] * (a ** 3 + a * a * b + a * b * b + b ** 3) / 4);
}

function insertInfluenceCut(cuts: number[], xi: number): void {
    if (xi < -INFLUENCE_ROOT_TOL || xi > 1 + INFLUENCE_ROOT_TOL) return;
    const clamped = xi < 0 ? 0 : xi > 1 ? 1 : xi;
    for (const existing of cuts) {
        if (Math.abs(existing - clamped) <= INFLUENCE_ROOT_TOL) return;
    }
    cuts.push(clamped);
}

function truckStationaryCuts(coeff: ArrayLike<number>): number[] {
    const cuts: number[] = [0, 1];
    let scale = 0;
    for (let i = 0; i < 4; i++) scale = Math.max(scale, Math.abs(coeff[i]));
    if (scale === 0) return cuts;
    const qa = 3 * coeff[3] / scale;
    const qb = 2 * coeff[2] / scale;
    const qc = coeff[1] / scale;
    if (Math.abs(qa) <= INFLUENCE_VALUE_TOL) {
        if (Math.abs(qb) > INFLUENCE_VALUE_TOL) insertInfluenceCut(cuts, -qc / qb);
    } else {
        let disc = qb * qb - 4 * qa * qc;
        if (Math.abs(disc) <= INFLUENCE_VALUE_TOL * (qb * qb + Math.abs(4 * qa * qc))) disc = 0;
        if (disc >= 0) {
            const q = qb >= 0 ? -0.5 * (qb + Math.sqrt(disc)) : -0.5 * (qb - Math.sqrt(disc));
            if (q === 0) {
                insertInfluenceCut(cuts, -qb / (2 * qa));
            } else {
                insertInfluenceCut(cuts, q / qa);
                insertInfluenceCut(cuts, qc / q);
            }
        }
    }
    return cuts;
}

function influenceCuts(coeff: ArrayLike<number>): number[] {
    let scale = 0;
    for (let k = 0; k < 4; k++) scale = Math.max(scale, Math.abs(coeff[k]));
    const cuts: number[] = [0, 1];
    if (scale === 0) return cuts;
    const stationary = truckStationaryCuts(coeff);
    stationary.sort((a, b) => a - b);
    for (const s of stationary) {
        if (Math.abs(influenceValue(coeff, s) / scale) <= INFLUENCE_VALUE_TOL) insertInfluenceCut(cuts, s);
    }
    for (let i = 0; i < stationary.length - 1; i++) {
        let a = stationary[i];
        let b = stationary[i + 1];
        let fa = influenceValue(coeff, a) / scale;
        const fb = influenceValue(coeff, b) / scale;
        if (Math.abs(fa) > INFLUENCE_VALUE_TOL && Math.abs(fb) > INFLUENCE_VALUE_TOL && fa * fb < 0) {
            for (let iter = 0; iter < 64; iter++) {
                const m = (a + b) / 2;
                const fm = influenceValue(coeff, m) / scale;
                if (fm === 0) { a = m; b = m; break; }
                else if (fa * fm < 0) b = m;
                else { a = m; fa = fm; }
                if (b - a <= INFLUENCE_ROOT_TOL) break;
            }
            insertInfluenceCut(cuts, (a + b) / 2);
        }
    }
    cuts.sort((a, b) => a - b);
    return cuts;
}

function truckLoadElement(position: number, side: number, nodes: ArrayLike<number>): number {
    const n = nodes.length - 1;
    const tol = INFLUENCE_ROOT_TOL * (1 + nodes[n]);
    if (position < -tol || position > nodes[n] + tol) return -1;
    if (Math.abs(position) <= tol && side < 0) return -1;
    if (Math.abs(position - nodes[n]) <= tol && side > 0) return -1;
    let lo = 0;
    let hi = n;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (nodes[mid] < position - tol) lo = mid + 1;
        else hi = mid;
    }
    let elem: number;
    if (Math.abs(nodes[lo] - position) <= tol && side > 0) elem = lo;
    else elem = lo - 1;
    if (elem < 0) elem = 0;
    if (elem >= n) elem = n - 1;
    return elem;
}

type TruckGeometry = {
    positions: number[];
    midElement: Int32Array;
    pointElement: Int32Array;
    z: Float64Array;
    delta: Float64Array;
    pointXi: Float64Array;
    intervals: number;
};

function uniqueSortedPositions(values: number[], tol: number): number[] {
    const sorted = [...values].sort((a, b) => a - b);
    const unique: number[] = [];
    for (const v of sorted) {
        if (unique.length === 0 || v - unique[unique.length - 1] > tol) unique.push(v);
    }
    return unique;
}

function buildTruckGeometry(nodes: number[], offsets: number[], nAxles: number): TruckGeometry {
    const raw: number[] = [];
    for (let p = 0; p < nodes.length; p++) {
        for (let i = 0; i < nAxles; i++) raw.push(nodes[p] + offsets[i]);
    }
    const tol = INFLUENCE_ROOT_TOL * (1 + nodes[nodes.length - 1]);
    const positions = uniqueSortedPositions(raw, tol);
    const intervals = Math.max(0, positions.length - 1);
    const midElement = new Int32Array(Math.max(0, nAxles * intervals));
    const pointElement = new Int32Array(Math.max(0, nAxles * intervals));
    const z = new Float64Array(Math.max(0, nAxles * intervals));
    const delta = new Float64Array(Math.max(0, nAxles * intervals));
    const pointXi = new Float64Array(Math.max(0, nAxles * intervals));
    if (positions.length < 2) return { positions, midElement, pointElement, z, delta, pointXi, intervals };
    for (let p = 0; p < positions.length - 1; p++) {
        const a = positions[p];
        const b = positions[p + 1];
        for (let i = 0; i < nAxles; i++) {
            const idx = i * intervals + p;
            const e = truckLoadElement((a + b) / 2 - offsets[i], 0, nodes);
            midElement[idx] = e;
            if (e >= 0) {
                z[idx] = (a - offsets[i] - nodes[e]) / (nodes[e + 1] - nodes[e]);
                delta[idx] = (b - a) / (nodes[e + 1] - nodes[e]);
            }
            const pe = truckLoadElement(a - offsets[i], 0, nodes);
            pointElement[idx] = pe;
            if (pe >= 0) {
                let xi = (a - offsets[i] - nodes[pe]) / (nodes[pe + 1] - nodes[pe]);
                if (xi < 0) xi = 0;
                if (xi > 1) xi = 1;
                pointXi[idx] = xi;
            }
        }
    }
    return { positions, midElement, pointElement, z, delta, pointXi, intervals };
}

function buildInfluenceCache(system: BeamSystem): { coeffs: Float64Array; nResponses: number } {
    const nElems = system.nElems;
    const nResponses = system.response.length;
    const coeffs = new Float64Array(nElems * nResponses * 4);
    const xs = [0.125, 0.375, 0.625, 0.875];
    for (let e = 0; e < nElems; e++) {
        const samples: Float64Array[] = [];
        for (let p = 0; p < 4; p++) {
            samples.push(system.unitResponse(e, xs[p]));
        }
        for (let r = 0; r < nResponses; r++) {
            const fitted = fitInfluenceCubic(samples[0][r], samples[1][r], samples[2][r], samples[3][r]);
            const base = (e * nResponses + r) * 4;
            coeffs[base] = fitted[0];
            coeffs[base + 1] = fitted[1];
            coeffs[base + 2] = fitted[2];
            coeffs[base + 3] = fitted[3];
        }
    }
    return { coeffs, nResponses };
}

function calculateUDLEnvelopes(
    coeffs: Float64Array, nElems: number, nResponses: number,
    elemLens: ArrayLike<number>, w: number
): { max: Float64Array; min: Float64Array } {
    const max = new Float64Array(nResponses);
    const min = new Float64Array(nResponses);
    if (w <= 0) return { max, min };
    const coeff = [0, 0, 0, 0];
    for (let e = 0; e < nElems; e++) {
        for (let r = 0; r < nResponses; r++) {
            const base = (e * nResponses + r) * 4;
            coeff[0] = coeffs[base];
            coeff[1] = coeffs[base + 1];
            coeff[2] = coeffs[base + 2];
            coeff[3] = coeffs[base + 3];
            const cuts = influenceCuts(coeff);
            for (let i = 0; i < cuts.length - 1; i++) {
                const area = w * elemLens[e] * influenceIntegral(coeff, cuts[i], cuts[i + 1]);
                if (area > 0) max[r] += area;
                else min[r] += area;
            }
        }
    }
    return { max, min };
}

function buildUdlTracer(
    system: BeamSystem, coeffs: Float64Array, nResponses: number, w: number,
    udlMax: Float64Array, udlMin: Float64Array
): UdlTracer {
    const middle = system.xNodes[system.nElems] / 2;
    let node = 0;
    for (let i = 1; i < system.nNodes; i++) {
        if (Math.abs(system.xNodes[i] - middle) < Math.abs(system.xNodes[node] - middle)) node = i;
    }
    const response = system.momentOffset + node;
    const intervals: UdlInterval[] = [];
    const influenceLine: UdlTracer['influenceLine'] = [];
    for (let e = 0; e < system.nElems; e++) {
        const base = (e * nResponses + response) * 4;
        const coeff = coeffs.subarray(base, base + 4);
        const cuts = influenceCuts(coeff);
        for (let i = 0; i < cuts.length - 1; i++) {
            const a = cuts[i];
            const b = cuts[i + 1];
            const contribution = w * system.lengths[e] * influenceIntegral(coeff, a, b);
            if (contribution !== 0) {
                intervals.push({
                    envelope: contribution > 0 ? 'max' : 'min', element: e + 1,
                    xiStart: a, xiEnd: b,
                    start: system.xNodes[e] + system.lengths[e] * a,
                    end: system.xNodes[e] + system.lengths[e] * b,
                    contribution,
                });
            }
        }
        const plotCuts = [...cuts];
        for (let i = 0; i <= 8; i++) insertInfluenceCut(plotCuts, i / 8);
        plotCuts.sort((a, b) => a - b);
        for (const xi of plotCuts) {
            const ordinate = influenceValue(coeff, xi);
            influenceLine.push({
                x: system.xNodes[e] + system.lengths[e] * xi, ordinate,
                max: Math.max(0, ordinate), min: Math.min(0, ordinate),
            });
        }
    }
    const max = udlMax[response];
    const min = udlMin[response];
    const reconstructedMax = system.partialUdlResponse(intervals, 'max', w)[response];
    const reconstructedMin = system.partialUdlResponse(intervals, 'min', w)[response];
    const tolerance = 1e-9 + UDL_VERIFY_TOL * Math.max(Math.abs(max), Math.abs(min));
    if (!Number.isFinite(max) || !Number.isFinite(min) ||
        Math.abs(reconstructedMax - max) > tolerance || Math.abs(reconstructedMin - min) > tolerance) {
        throw new Error(`Partial-UDL FEM reconstruction failed at moment node ${node + 1}: ` +
            `max ${reconstructedMax} versus ${max}; min ${reconstructedMin} versus ${min}.`);
    }
    return {
        response: 'moment', nodeIndex: node + 1, x: system.xNodes[node], unit: 'kNm',
        intensity: w, max, min, reconstructedMax, reconstructedMin, tolerance,
        intervals, influenceLine,
    };
}

type PolyResult = { poly: Float64Array; masks: Int32Array; bases: Float64Array };

function truckPolynomialPair(
    response: number, a: number, b: number, side: number,
    weights: ArrayLike<number>, offsets: ArrayLike<number>, nAxles: number,
    nodes: ArrayLike<number>, coeffs: Float64Array, nResponses: number, nElems: number,
    autoDLA: boolean, multiplier: number, axleFactor: number, isReverse: boolean,
    cursor: Int32Array, responseCoeffs: Float64Array | null, useCursor: boolean,
    useCoefficients: boolean, includePoint: boolean, geometry: TruckGeometry | null,
    geometryIndex: number
): PolyResult {
    const poly = new Float64Array(16);
    const masks = new Int32Array(4);
    const bases = new Float64Array(4);
    const counts = [0, 0, 0, 0];
    let p0 = 0, p1 = 0, p2 = 0, p3 = 0;
    let q0 = 0, q1 = 0, q2 = 0, q3 = 0;
    let pointMax = 0, pointMin = 0;
    const n = nodes.length - 1;
    const totalLength = nodes[n];
    const tol = INFLUENCE_ROOT_TOL * (1 + totalLength);
    const lastRow = includePoint ? 3 : 1;
    let bit = isReverse ? Math.pow(2, nAxles - 1) : 1;
    const intervals = geometry?.intervals ?? 0;
    for (let i = 0; i < nAxles; i++) {
        const position = a - offsets[i];
        let e = -1;
        if (geometry && geometryIndex >= 0) {
            e = geometry.midElement[i * intervals + geometryIndex];
            if (useCursor) cursor[i] = e;
        } else if (includePoint) {
            const midpoint = (a + b) / 2 - offsets[i];
            e = -1;
            if (midpoint >= -tol && midpoint <= totalLength + tol) {
                e = cursor[i];
                if (e < 0) e = 0;
                const left = midpoint - tol;
                while (e < nElems - 1 && nodes[e + 1] < left) e++;
                while (e > 0 && nodes[e] >= left) e--;
                // Clamp into range
                if (e < 0) e = 0;
                if (e >= nElems) e = nElems - 1;
            }
            cursor[i] = e;
        } else if (useCursor) {
            e = truckLoadElementWithHint((a + b) / 2 - offsets[i], side, nodes, cursor[i]);
            cursor[i] = e;
        } else {
            e = truckLoadElement((a + b) / 2 - offsets[i], side, nodes);
        }
        if (e >= 0 && weights[i] !== 0) {
            let curZ: number;
            let curDelta: number;
            if (geometry && geometryIndex >= 0) {
                curZ = geometry.z[i * intervals + geometryIndex];
                curDelta = geometry.delta[i * intervals + geometryIndex];
            } else {
                curZ = (a - offsets[i] - nodes[e]) / (nodes[e + 1] - nodes[e]);
                curDelta = (b - a) / (nodes[e + 1] - nodes[e]);
            }
            let zz = curZ;
            if (a === b) {
                if (zz < 0) zz = 0;
                if (zz > 1) zz = 1;
            }
            let c0: number, c1: number, c2: number, c3: number, scale: number;
            if (useCoefficients && responseCoeffs) {
                c0 = responseCoeffs[e * 7];
                c1 = responseCoeffs[e * 7 + 1];
                c2 = responseCoeffs[e * 7 + 2];
                c3 = responseCoeffs[e * 7 + 3];
                scale = responseCoeffs[e * 7 + 4];
            } else {
                const base = (e * nResponses + response) * 4;
                c0 = coeffs[base];
                c1 = coeffs[base + 1];
                c2 = coeffs[base + 2];
                c3 = coeffs[base + 3];
                scale = Math.abs(c0);
                if (Math.abs(c1) > scale) scale = Math.abs(c1);
                if (Math.abs(c2) > scale) scale = Math.abs(c2);
                if (Math.abs(c3) > scale) scale = Math.abs(c3);
            }
            const xiMid = zz + curDelta / 2;
            const unitMid = ((c3 * xiMid + c2) * xiMid + c1) * xiMid + c0;
            const valueMid = weights[i] * unitMid;
            const threshold = INFLUENCE_VALUE_TOL * scale * Math.abs(weights[i]);
            let row = -1;
            if (valueMid > threshold) row = 0;
            if (valueMid < -threshold) row = 1;
            if (row >= 0) {
                if (a === b) {
                    if (row === 0) p0 += valueMid;
                    else q0 += valueMid;
                } else {
                    const unit0 = ((c3 * zz + c2) * zz + c1) * zz + c0;
                    if (row === 0) {
                        p0 += weights[i] * unit0;
                        p1 += weights[i] * curDelta * (c1 + 2 * c2 * zz + 3 * c3 * zz * zz);
                        p2 += weights[i] * curDelta * curDelta * (c2 + 3 * c3 * zz);
                        p3 += weights[i] * Math.pow(curDelta, 3) * c3;
                    } else {
                        q0 += weights[i] * unit0;
                        q1 += weights[i] * curDelta * (c1 + 2 * c2 * zz + 3 * c3 * zz * zz);
                        q2 += weights[i] * curDelta * curDelta * (c2 + 3 * c3 * zz);
                        q3 += weights[i] * Math.pow(curDelta, 3) * c3;
                    }
                }
                masks[row] |= bit;
                counts[row]++;
            }
        }
        if (includePoint && weights[i] !== 0) {
            let pe = -1;
            if (geometry && geometryIndex >= 0) {
                pe = geometry.pointElement[i * intervals + geometryIndex];
            } else {
                pe = -1;
                if (position >= -tol && position <= totalLength + tol) {
                    pe = e;
                    if (pe < 0) pe = nElems - 1;
                    const left = position - tol;
                    while (pe > 0 && nodes[pe] >= left) pe--;
                }
            }
            if (pe >= 0) {
                let pc0: number, pc1: number, pc2: number, pc3: number, pscale: number;
                if (pe !== e || (geometry && geometryIndex >= 0)) {
                    if (useCoefficients && responseCoeffs) {
                        pc0 = responseCoeffs[pe * 7];
                        pc1 = responseCoeffs[pe * 7 + 1];
                        pc2 = responseCoeffs[pe * 7 + 2];
                        pc3 = responseCoeffs[pe * 7 + 3];
                        pscale = responseCoeffs[pe * 7 + 4];
                    } else {
                        const base = (pe * nResponses + response) * 4;
                        pc0 = coeffs[base];
                        pc1 = coeffs[base + 1];
                        pc2 = coeffs[base + 2];
                        pc3 = coeffs[base + 3];
                        pscale = Math.abs(pc0);
                        if (Math.abs(pc1) > pscale) pscale = Math.abs(pc1);
                        if (Math.abs(pc2) > pscale) pscale = Math.abs(pc2);
                        if (Math.abs(pc3) > pscale) pscale = Math.abs(pc3);
                    }
                } else {
                    // Reuse current element coefficients
                    if (useCoefficients && responseCoeffs) {
                        pc0 = responseCoeffs[e * 7];
                        pc1 = responseCoeffs[e * 7 + 1];
                        pc2 = responseCoeffs[e * 7 + 2];
                        pc3 = responseCoeffs[e * 7 + 3];
                        pscale = responseCoeffs[e * 7 + 4];
                    } else {
                        const base = (e * nResponses + response) * 4;
                        pc0 = coeffs[base];
                        pc1 = coeffs[base + 1];
                        pc2 = coeffs[base + 2];
                        pc3 = coeffs[base + 3];
                        pscale = Math.abs(pc0);
                        if (Math.abs(pc1) > pscale) pscale = Math.abs(pc1);
                        if (Math.abs(pc2) > pscale) pscale = Math.abs(pc2);
                        if (Math.abs(pc3) > pscale) pscale = Math.abs(pc3);
                    }
                }
                let pxi: number;
                if (geometry && geometryIndex >= 0) pxi = geometry.pointXi[i * intervals + geometryIndex];
                else {
                    pxi = (a - offsets[i] - nodes[pe]) / (nodes[pe + 1] - nodes[pe]);
                    if (pxi < 0) pxi = 0;
                    if (pxi > 1) pxi = 1;
                }
                const unitP = ((pc3 * pxi + pc2) * pxi + pc1) * pxi + pc0;
                const valueP = weights[i] * unitP;
                const thresholdP = INFLUENCE_VALUE_TOL * pscale * Math.abs(weights[i]);
                let prow = -1;
                if (valueP > thresholdP) prow = 2;
                if (valueP < -thresholdP) prow = 3;
                if (prow >= 0) {
                    if (prow === 2) pointMax += valueP;
                    else pointMin += valueP;
                    masks[prow] |= bit;
                    counts[prow]++;
                }
            }
        }
        if (isReverse) bit = Math.floor(bit / 2);
        else bit = bit * 2;
    }
    for (let row = 0; row <= lastRow; row++) {
        bases[row] = 0;
        if (autoDLA) bases[row] = truckGroupDla(counts[row], masks[row] === 7);
        const factor = axleFactor * (1 + bases[row] * multiplier);
        let v: number;
        if (row === 0) v = p0;
        else if (row === 1) v = q0;
        else if (row === 2) v = pointMax;
        else v = pointMin;
        if (a === b || row >= 2) {
            poly[row * 4] = v * axleFactor * (1 + bases[row] * multiplier);
            poly[row * 4 + 1] = 0;
            poly[row * 4 + 2] = 0;
            poly[row * 4 + 3] = 0;
        } else {
            poly[row * 4] = v * factor;
            if (row === 0) {
                poly[row * 4 + 1] = p1 * factor;
                poly[row * 4 + 2] = p2 * factor;
                poly[row * 4 + 3] = p3 * factor;
            } else {
                poly[row * 4 + 1] = q1 * factor;
                poly[row * 4 + 2] = q2 * factor;
                poly[row * 4 + 3] = q3 * factor;
            }
        }
    }
    return { poly, masks, bases };
}

function truckLoadElementWithHint(position: number, side: number, nodes: ArrayLike<number>, hint: number): number {
    const n = nodes.length - 1;
    const tol = INFLUENCE_ROOT_TOL * (1 + nodes[n]);
    if (position < -tol || position > nodes[n] + tol) return -1;
    if (Math.abs(position) <= tol && side < 0) return -1;
    if (Math.abs(position - nodes[n]) <= tol && side > 0) return -1;
    let lo = 0;
    let hi = n;
    if (hint >= 0 && hint < n && nodes[hint] < position - tol) {
        lo = hint;
        while (lo < n && nodes[lo] < position - tol) lo++;
        hi = lo;
    } else {
        while (lo < hi) {
            const mid = (lo + hi) >> 1;
            if (nodes[mid] < position - tol) lo = mid + 1;
            else hi = mid;
        }
    }
    let elem: number;
    if (Math.abs(nodes[lo] - position) <= tol && side > 0) elem = lo;
    else elem = lo - 1;
    if (elem < 0) elem = 0;
    if (elem >= n) elem = n - 1;
    return elem;
}

function prepareTruckResponse(
    response: number, nodes: number[], coeffs: Float64Array,
    nElems: number, nResponses: number
): { knots: number[]; responseCoeffs: Float64Array; isZero: boolean } {
    const knots: number[] = [];
    const responseCoeffs = new Float64Array(nElems * 7);
    let isZero = true;
    const totalLength = nodes[nodes.length - 1];
    const coeff = [0, 0, 0, 0];
    for (let e = 0; e < nElems; e++) {
        const base = (e * nResponses + response) * 4;
        let scale = 0;
        for (let p = 0; p < 4; p++) {
            coeff[p] = coeffs[base + p];
            responseCoeffs[e * 7 + p] = coeff[p];
            scale = Math.max(scale, Math.abs(coeff[p]));
        }
        responseCoeffs[e * 7 + 4] = scale;
        if (scale !== 0) isZero = false;
        let upper = coeff[0];
        let lower = coeff[0];
        for (let p = 1; p < 4; p++) {
            if (coeff[p] > 0) upper += coeff[p];
            else lower += coeff[p];
        }
        const extension = 4 * INFLUENCE_ROOT_TOL * (1 + totalLength) / (nodes[e + 1] - nodes[e]);
        const allowance = extension * (Math.abs(coeff[1]) + (2 + extension) * Math.abs(coeff[2]) +
            (3 + 3 * extension + extension * extension) * Math.abs(coeff[3])) +
            INFLUENCE_VALUE_TOL * (1 + 4 * scale);
        responseCoeffs[e * 7 + 5] = upper + allowance;
        responseCoeffs[e * 7 + 6] = lower - allowance;
        const cuts = influenceCuts(coeff);
        for (const c of cuts) {
            const pos = nodes[e] + c * (nodes[e + 1] - nodes[e]);
            if (knots.length === 0 || pos !== knots[knots.length - 1]) knots.push(pos);
        }
    }
    return { knots, responseCoeffs, isZero };
}

function truckGeometryCanImprove(
    index: number, geometry: TruckGeometry, responseCoeffs: Float64Array,
    weights: ArrayLike<number>, nAxles: number, factor: number, bestMax: number, bestMin: number
): boolean {
    let maxValue = 0;
    let minValue = 0;
    const intervals = geometry.intervals;
    for (let i = 0; i < nAxles; i++) {
        const e = geometry.midElement[i * intervals + index];
        const point = geometry.pointElement[i * intervals + index];
        let upper = 0;
        let lower = 0;
        if (e >= 0) {
            const up = responseCoeffs[e * 7 + 5];
            const lo = responseCoeffs[e * 7 + 6];
            if (up > upper) upper = up;
            if (lo < lower) lower = lo;
        }
        if (point >= 0 && point !== e) {
            const up = responseCoeffs[point * 7 + 5];
            const lo = responseCoeffs[point * 7 + 6];
            if (up > upper) upper = up;
            if (lo < lower) lower = lo;
        }
        const w = weights[i];
        if (w >= 0) {
            maxValue += w * upper;
            minValue += w * lower;
        } else {
            maxValue += w * lower;
            minValue += w * upper;
        }
    }
    maxValue *= factor;
    minValue *= factor;
    const margin = INFLUENCE_VALUE_TOL * (1 + Math.abs(maxValue) + Math.abs(minValue));
    return maxValue + margin > bestMax || minValue - margin < bestMin;
}

function truckResponseExtrema(
    response: number, weights: number[], offsets: number[], nAxles: number,
    nodes: number[], knots: number[], responseCoeffs: Float64Array,
    geometry: TruckGeometry, autoDLA: boolean, multiplier: number,
    axleFactor: number, isReverse: boolean,
    best: number[], govLead: number[], govMask: Int32Array,
    govBase: number[], govSide: Int32Array,
    coeffs: Float64Array, nResponses: number, nElems: number
): void {
    const raw: number[] = [];
    for (const k of knots) {
        for (let i = 0; i < nAxles; i++) raw.push(k + offsets[i]);
    }
    const tol = INFLUENCE_ROOT_TOL * (1 + nodes[nodes.length - 1]);
    const positions = uniqueSortedPositions(raw, tol);
    const unique = positions.length;
    if (unique === 0) return;
    const geometryLast = geometry.positions.length - 1;
    let boundFactor = axleFactor;
    if (autoDLA) boundFactor = boundFactor * (1 + 0.4 * multiplier);
    const cursor = new Int32Array(nAxles);
    let geometryPosition = 0;
    for (let i = 0; i < unique; i++) {
        const a = positions[i];
        let b = a;
        let pointOffset = 0;
        if (i < unique - 1) {
            b = positions[i + 1];
            pointOffset = 2;
        }
        while (geometryPosition < geometryLast && geometry.positions[geometryPosition] < a) geometryPosition++;
        let geometryIndex = -1;
        if (pointOffset === 2 && geometryPosition < geometryLast) {
            if (a === geometry.positions[geometryPosition] && b === geometry.positions[geometryPosition + 1]) {
                geometryIndex = geometryPosition;
            }
        }
        if (geometryIndex >= 0 && boundFactor >= 0) {
            if (!truckGeometryCanImprove(geometryIndex, geometry, responseCoeffs, weights, nAxles, boundFactor, best[0], best[1])) {
                continue;
            }
        }
        const { poly, masks, bases } = truckPolynomialPair(
            response, a, b, 0, weights, offsets, nAxles, nodes, coeffs, nResponses, nElems,
            autoDLA, multiplier, axleFactor, isReverse, cursor, responseCoeffs, true, true,
            pointOffset === 2, geometry, geometryIndex
        );
        for (let row = 0; row <= 1; row++) {
            const value = poly[(row + pointOffset) * 4];
            if ((row === 0 && value > best[row]) || (row === 1 && value < best[row])) {
                best[row] = value;
                govLead[row] = a;
                govMask[row] = masks[row + pointOffset];
                govBase[row] = bases[row + pointOffset];
                govSide[row] = 0;
            }
        }
        if (i < unique - 1) {
            for (let row = 0; row <= 1; row++) {
                if (masks[row] === 0) continue;
                const c0 = poly[row * 4];
                const c1 = poly[row * 4 + 1];
                const c2 = poly[row * 4 + 2];
                const c3 = poly[row * 4 + 3];
                let bound = c0;
                let boundScale = Math.abs(c0);
                const powers = [c1, c2, c3];
                for (const cp of powers) {
                    boundScale += Math.abs(cp);
                    if ((row === 0 && cp > 0) || (row === 1 && cp < 0)) bound += cp;
                }
                const margin = INFLUENCE_VALUE_TOL * (1 + boundScale);
                if ((row === 0 && bound + margin <= best[row]) ||
                    (row === 1 && bound - margin >= best[row])) continue;
                const d0 = c1;
                const d1 = c1 + c2;
                const d2 = c1 + 2 * c2 + 3 * c3;
                let candidates = [0, 1];
                if (!((d0 >= 0 && d1 >= 0 && d2 >= 0) || (d0 <= 0 && d1 <= 0 && d2 <= 0))) {
                    const stationary = truckStationaryCuts([c0, c1, c2, c3]);
                    // truckStationaryCuts returns [0,1,...]; use interior points
                    candidates = stationary;
                }
                for (const xi of candidates) {
                    const lead = a + xi * (b - a);
                    const value = ((c3 * xi + c2) * xi + c1) * xi + c0;
                    if ((row === 0 && value > best[row]) || (row === 1 && value < best[row])) {
                        let side = 0;
                        if (xi === 0) side = 1;
                        if (xi === 1) side = -1;
                        best[row] = value;
                        govLead[row] = lead;
                        govMask[row] = masks[row];
                        govBase[row] = bases[row];
                        govSide[row] = side;
                    }
                }
            }
        }
    }
}

function runInfluenceTruckEnvelope(
    axleFactor: number, autoDLA: boolean, multiplier: number,
    baseAxles: number[], spacings: number[], nAxles: number,
    nodes: number[], nElems: number, nResponses: number,
    coeffs: Float64Array, steps: number[], nSupports: number,
    udlMax: Float64Array, udlMin: Float64Array,
    onProgress?: (fraction: number, message: string) => void,
    caseLabel = 'Truck'
): {
    maxima: Float64Array; minima: Float64Array;
    histMax: Float64Array; histMin: Float64Array;
    optMax: Float64Array; optMin: Float64Array; optGov: Float64Array;
    truckSolves: number;
} {
    const maxima = new Float64Array(nResponses).fill(0);
    const minima = new Float64Array(nResponses).fill(0);
    const rOffset = 4 * nElems + 2;
    // Forward / reverse weights and offsets
    const forwardWeights = [...baseAxles];
    const forwardOffsets = new Array(nAxles).fill(0);
    for (let i = 1; i < nAxles; i++) forwardOffsets[i] = forwardOffsets[i - 1] + spacings[i - 1];
    const reverseWeights = [...baseAxles].reverse();
    const reverseOffsets = new Array(nAxles).fill(0);
    for (let i = 1; i < nAxles; i++) reverseOffsets[i] = reverseOffsets[i - 1] + spacings[nAxles - 1 - i];
    const forwardGeometry = buildTruckGeometry(nodes, forwardOffsets, nAxles);
    const reverseGeometry = buildTruckGeometry(nodes, reverseOffsets, nAxles);
    const optGov = new Float64Array(nSupports).fill(0);
    let truckSolves = 0;
    const best = [0, 0];
    const govLead = [0, 0];
    const govMask = new Int32Array(2);
    const govBase = [0, 0];
    const govSide = new Int32Array(2);
    for (let r = 0; r < nResponses; r++) {
        if (r % 20 === 0) onProgress?.(r / nResponses * 0.85, `${caseLabel}: optimising truck response ${r + 1} of ${nResponses}`);
        const { knots, responseCoeffs, isZero } = prepareTruckResponse(r, nodes, coeffs, nElems, nResponses);
        if (isZero) continue;
        best[0] = maxima[r];
        best[1] = minima[r];
        const oldMax = best[0];
        // Forward
        truckResponseExtrema(r, forwardWeights, forwardOffsets, nAxles, nodes, knots, responseCoeffs,
            forwardGeometry, autoDLA, multiplier, axleFactor, false, best, govLead, govMask, govBase, govSide,
            coeffs, nResponses, nElems);
        // Reverse
        truckResponseExtrema(r, reverseWeights, reverseOffsets, nAxles, nodes, knots, responseCoeffs,
            reverseGeometry, autoDLA, multiplier, axleFactor, true, best, govLead, govMask, govBase, govSide,
            coeffs, nResponses, nElems);
        maxima[r] = best[0];
        minima[r] = best[1];
        if (r >= rOffset && best[0] > oldMax) optGov[r - rOffset] = govLead[0];
        truckSolves += 2;
    }
    // Sampled reaction histories at sweep steps (for diagrams), with selected-group DLA
    const nSteps = steps.length;
    const histMax = new Float64Array(nSteps * nSupports).fill(-Infinity);
    const histMin = new Float64Array(nSteps * nSupports).fill(Infinity);
    const emptyGeometry: TruckGeometry = {
        positions: [], midElement: new Int32Array(0), pointElement: new Int32Array(0),
        z: new Float64Array(0), delta: new Float64Array(0), pointXi: new Float64Array(0), intervals: 0,
    };
    for (let direction = 0; direction < 2; direction++) {
        const weights = direction === 0 ? forwardWeights : reverseWeights;
        const offsets = direction === 0 ? forwardOffsets : reverseOffsets;
        const isReverse = direction === 1;
        const cursor = new Int32Array(nAxles);
        for (let p = 0; p < nSteps; p++) {
            for (let s = 0; s < nSupports; s++) {
                const r = rOffset + s;
                const { poly } = truckPolynomialPair(
                    r, steps[p], steps[p], 0, weights, offsets, nAxles, nodes, coeffs,
                    nResponses, nElems, autoDLA, multiplier, axleFactor, isReverse,
                    cursor, null, true, false, false, emptyGeometry, -1
                );
                const hi = poly[0] + udlMax[rOffset + s];
                const lo = poly[1 * 4] + udlMin[rOffset + s];
                const idx = p * nSupports + s;
                if (hi > histMax[idx]) histMax[idx] = hi;
                if (lo < histMin[idx]) histMin[idx] = lo;
            }
            truckSolves++;
            if (p % 50 === 0) onProgress?.(0.85 + 0.15 * (direction * nSteps + p) / (nSteps * 2), `${caseLabel}: reaction curves ${p + 1}/${nSteps}`);
        }
    }
    // Combine continuous optima with UDL envelopes
    const optMax = new Float64Array(nSupports);
    const optMin = new Float64Array(nSupports);
    for (let s = 0; s < nSupports; s++) {
        optMax[s] = maxima[rOffset + s] + udlMax[rOffset + s];
        optMin[s] = minima[rOffset + s] + udlMin[rOffset + s];
    }
    // Add UDL to V/M/D optima (histories already include UDL)
    // maxima/minima arrays are truck-only; final envelopes add UDL at packaging time.
    return { maxima, minima, histMax, histMin, optMax, optMin, optGov, truckSolves };
}

export function analyzeBeam(
    request: AnalysisRequest, onProgress?: (progress: AnalysisProgress) => void
): AnalysisResults {
    const started = performance.now();
    validateInputs(request);
    const { spans, axles, config } = request;
    const dla = resolveDla(config);
    onProgress?.({ fraction: 0, message: 'Factorizing beam stiffness...' });
    const system = new BeamSystem(spans, config);
    const nElems = system.nElems;
    const nSupports = system.supports.length;
    const nResponses = system.response.length;
    const stepInfo = computeEffectiveIncrement(spans, axles, config.truckIncrement, config.nElemsPerSpan);
    const positions = buildTruckPositions(system.supports, axles, stepInfo.effective);
    if (positions.length * nSupports > 2000000)
        throw new Error('Reaction histories exceed 2,000,000 ordinates per case. Reduce the span or axle count.');
    onProgress?.({ fraction: 0.02, message: 'Generating influence functions...' });
    const { coeffs } = buildInfluenceCache(system);
    let udlIntegrations = 0;
    const baseAxles = axles.map(a => a.load);
    const spacings = axles.map(a => a.spacing);
    const nAxles = axles.length;

    const runCase = (
        label: 'Truck' | 'Lane', axleFactor: number, wUdl: number,
        autoDLA: boolean, multiplier: number
    ) => {
        const { max: udlMax, min: udlMin } = calculateUDLEnvelopes(coeffs, nElems, nResponses, system.lengths, wUdl);
        if (wUdl > 0) udlIntegrations += nElems * nResponses;
        const env = runInfluenceTruckEnvelope(
            axleFactor, autoDLA, multiplier, baseAxles, spacings, nAxles,
            system.xNodes, nElems, nResponses, coeffs, positions, nSupports,
            udlMax, udlMin,
            (fraction, message) => onProgress?.({
                fraction: (config.loadCase === 'envelope' && label === 'Lane' ? 0.48 : 0.02) +
                    fraction * (config.loadCase === 'envelope' ? 0.46 : 0.92),
                message,
            }),
            label
        );
        // Final envelopes: continuous truck optima + UDL zones
        const finalMax = new Float64Array(nResponses);
        const finalMin = new Float64Array(nResponses);
        for (let i = 0; i < nResponses; i++) {
            finalMax[i] = env.maxima[i] + udlMax[i];
            finalMin[i] = env.minima[i] + udlMin[i];
        }
        return { ...env, udlMax, udlMin, finalMax, finalMin };
    };

    type BuiltCase = ReturnType<typeof runCase> & { dlaAuto: boolean; dlaBase: number; dlaMultiplier: number; dlaUsed: number };
    const built: Partial<Record<LoadCase, BuiltCase>> = {};
    if (config.loadCase !== 'lane') {
        const axleFactor = 1 + (dla.isAuto ? 0 : dla.effective);
        const t = runCase('Truck', axleFactor, 0, dla.isAuto, dla.multiplier);
        built.truck = { ...t, dlaAuto: dla.isAuto, dlaBase: dla.isAuto ? 0 : dla.base, dlaMultiplier: dla.multiplier, dlaUsed: dla.isAuto ? 0 : dla.effective };
    }
    if (config.loadCase !== 'truck') {
        const wLane = config.laneUdl ?? LANE_UDL;
        const l = runCase('Lane', LANE_TRUCK_FACTOR, wLane, false, 1);
        built.lane = { ...l, dlaAuto: false, dlaBase: 0, dlaMultiplier: 1, dlaUsed: 0 };
    }

    const makeCase = (b: BuiltCase): CaseResults => {
        const shear: EnvelopePoint[] = system.xShear.map((x, i) => ({ x, max: b.finalMax[i], min: b.finalMin[i] }));
        const moment: EnvelopePoint[] = system.xNodes.map((x, i) => ({ x, max: b.finalMax[2 * nElems + i], min: b.finalMin[2 * nElems + i] }));
        const deflection: EnvelopePoint[] = system.xNodes.map((x, i) => ({ x, max: b.finalMax[3 * nElems + 1 + i], min: b.finalMin[3 * nElems + 1 + i] }));
        const reactionDiagrams = system.supports.map((_, s) =>
            positions.map((x, p) => ({
                x, max: b.histMax[p * nSupports + s], min: b.histMin[p * nSupports + s],
            })));
        const reactions: ReactionEnvelope[] = system.supports.map((x, s) => ({
            x, max: b.optMax[s], min: b.optMin[s], govPos: b.optGov[s],
        }));
        return {
            shear, moment, deflection, reactionDiagrams, reactions,
            dlaUsed: b.dlaUsed, dlaAuto: b.dlaAuto, dlaBase: b.dlaBase, dlaMultiplier: b.dlaMultiplier,
        };
    };

    const cases: Partial<Record<LoadCase, CaseResults>> = {};
    if (built.truck) cases.truck = makeCase(built.truck);
    if (built.lane) cases.lane = makeCase(built.lane);
    if (cases.truck && cases.lane && built.truck && built.lane) {
        const truck = cases.truck;
        const lane = cases.lane;
        const bt = built.truck;
        const bl = built.lane;
        const combine = (t: EnvelopePoint[], l: EnvelopePoint[]) =>
            t.map((p, i) => ({ x: p.x, max: Math.max(p.max, l[i].max), min: Math.min(p.min, l[i].min) }));
        // Reaction histories: envelope of both cases per ordinate (VBA parity)
        const reactionDiagrams = truck.reactionDiagrams.map((diagram, s) =>
            diagram.map((p, idx) => ({
                x: p.x,
                max: Math.max(p.max, lane.reactionDiagrams[s][idx].max),
                min: Math.min(p.min, lane.reactionDiagrams[s][idx].min),
            })));
        // Support summary: optimised maxima envelope (VBA: max of optimised, min of min)
        const reactions = truck.reactions.map((p, s) => {
            const lmax = lane.reactions[s].max;
            const lmin = lane.reactions[s].min;
            const useLane = lmax > p.max;
            return {
                x: p.x,
                max: Math.max(p.max, lmax),
                min: Math.min(p.min, lmin),
                govPos: useLane ? lane.reactions[s].govPos : p.govPos,
            };
        });
        cases.envelope = {
            shear: combine(truck.shear, lane.shear),
            moment: combine(truck.moment, lane.moment),
            deflection: combine(truck.deflection, lane.deflection),
            reactionDiagrams,
            reactions,
            dlaUsed: bt.dlaUsed,
            dlaAuto: bt.dlaAuto,
            dlaBase: bt.dlaBase,
            dlaMultiplier: bt.dlaMultiplier,
        };
        void bl;
    }
    const selected = cases[config.loadCase];
    if (!selected) throw new Error('The selected load case was not calculated.');
    if (built.lane) onProgress?.({ fraction: 0.95, message: 'Verifying partial-UDL influence zones...' });
    const udlTracer = built.lane
        ? buildUdlTracer(system, coeffs, nResponses, config.laneUdl ?? LANE_UDL, built.lane.udlMax, built.lane.udlMin)
        : undefined;
    onProgress?.({ fraction: 1, message: 'Analysis complete.' });
    let truckSolves = 0;
    if (built.truck) truckSolves += built.truck.truckSolves;
    if (built.lane) truckSolves += built.lane.truckSolves;
    return {
        ...selected, cases, loadCase: config.loadCase, spans: spans.map(s => ({ ...s })),
        axles: axles.map(a => ({ ...a })), config: { ...config, dlaMultiplier: dla.multiplier, laneUdl: config.laneUdl ?? LANE_UDL },
        xNodes: system.xNodes, supportPositions: system.supports, truckPositions: positions,
        incrementUsed: stepInfo.effective, baseIncrement: config.truckIncrement,
        incrementReason: stepInfo.reason, elapsedMs: performance.now() - started,
        udlTracer,
        stats: { factorizations: 1, truckSolves, influenceSolves: 4 * nElems, udlIntegrations },
    };
}
