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
};
export type EnvelopePoint = { x: number; max: number; min: number };
export type ReactionEnvelope = EnvelopePoint & { govPos: number };
export type CaseResults = {
    shear: EnvelopePoint[];
    moment: EnvelopePoint[];
    deflection: EnvelopePoint[];
    reactions: ReactionEnvelope[];
    reactionDiagrams: EnvelopePoint[][];
    dlaUsed: number;
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
    stats: { factorizations: number; truckSolves: number; udlSolves: number };
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
    loadCase: 'truck', dlaOverride: null,
};
export const MAX_PATTERN_SPANS = 12;
export const MAX_AXLES = 20;
export const MAX_SWEEP_STEPS = 6000;
const MIN_SWEEP_STEP = 0.02;
const BAND = 3;

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
    const totalLength = spans.reduce((sum, s) => sum + s.length, 0);
    const truckLength = axles.slice(0, -1).reduce((sum, a) => sum + a.spacing, 0);
    if (!Number.isFinite(totalLength + 2 * truckLength) || !Number.isFinite(config.E * config.I))
        throw new Error('Geometry or stiffness exceeds the supported numeric range.');
    if (!Number.isSafeInteger(spans.length * config.nElemsPerSpan) ||
        spans.length * config.nElemsPerSpan > 10000)
        throw new Error('The model is limited to 10,000 beam elements. Reduce the mesh resolution.');
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
    private readonly lengths: Float64Array;
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
        const x2 = xi * xi;
        const x3 = x2 * xi;
        const shapes = [1 - 3 * x2 + 2 * x3, le * (xi - 2 * x2 + x3),
            3 * x2 - 2 * x3, le * (-x2 + x3)];
        for (let i = 0; i < 4; i++) {
            const value = -magnitude * shapes[i] * 1000;
            this.load[low * 2 + i] += value;
            this.elemLoads[low * 4 + i] += value;
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

    udl(spanIndex: number, intensity: number): Float64Array {
        this.clearLoads();
        const first = spanIndex * this.config.nElemsPerSpan;
        for (let e = first; e < first + this.config.nElemsPerSpan; e++) {
            const le = this.lengths[e];
            const force = -intensity * le / 2 * 1000;
            const moment = -intensity * le * le / 12 * 1000;
            const values = [force, moment, force, -moment];
            for (let i = 0; i < 4; i++) {
                this.load[e * 2 + i] += values[i];
                this.elemLoads[e * 4 + i] = values[i];
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

export function analyzeBeam(
    request: AnalysisRequest, onProgress?: (progress: AnalysisProgress) => void
): AnalysisResults {
    const started = performance.now();
    validateInputs(request);
    const { spans, axles, config } = request;
    onProgress?.({ fraction: 0, message: 'Factorizing beam stiffness...' });
    const system = new BeamSystem(spans, config);
    const stepInfo = computeEffectiveIncrement(spans, axles, config.truckIncrement, config.nElemsPerSpan);
    const positions = buildTruckPositions(system.supports, axles, stepInfo.effective);
    if (positions.length * system.supports.length > 2000000)
        throw new Error('Reaction histories exceed 2,000,000 ordinates per case. Reduce the span or axle count.');
    const max = new Float64Array(system.response.length).fill(-Infinity);
    const min = new Float64Array(system.response.length).fill(Infinity);
    const udlMax = new Float64Array(system.response.length);
    const udlMin = new Float64Array(system.response.length);
    let udlSolves = 0;
    if (config.loadCase !== 'truck') {
        // Every load pattern is a sum of independent span responses. Choosing each
        // positive/negative contribution gives exactly the same extrema as all 2^n patterns.
        for (let s = 0; s < spans.length; s++) {
            const response = system.udl(s, 9);
            for (let i = 0; i < response.length; i++) {
                udlMax[i] += Math.max(0, response[i]);
                udlMin[i] += Math.min(0, response[i]);
            }
            udlSolves++;
        }
    }
    const reactionMax = positions.map(() => new Float64Array(system.supports.length).fill(-Infinity));
    const reactionMin = positions.map(() => new Float64Array(system.supports.length).fill(Infinity));
    let truckSolves = 0;
    for (const [directionIndex, direction] of [axles, reverseAxles(axles)].entries()) {
        for (let p = 0; p < positions.length; p++) {
            const response = system.truck(positions[p], direction);
            for (let i = 0; i < response.length; i++) {
                max[i] = Math.max(max[i], response[i]);
                min[i] = Math.min(min[i], response[i]);
            }
            for (let s = 0; s < system.supports.length; s++) {
                const value = response[system.reactionOffset + s];
                reactionMax[p][s] = Math.max(reactionMax[p][s], value);
                reactionMin[p][s] = Math.min(reactionMin[p][s], value);
            }
            truckSolves++;
            if (p % 25 === 0) onProgress?.({
                fraction: (directionIndex * positions.length + p) / (positions.length * 2),
                message: `${directionIndex === 0 ? 'Forward' : 'Reverse'} truck sweep: ${positions[p].toFixed(2)}m`,
            });
        }
    }
    const dla = config.dlaOverride ?? computeAutoDlaInfo(spans, axles).dla;
    const makeCase = (factor: number, lane: boolean): CaseResults => {
        const hi = (i: number) => max[i] * factor + (lane ? udlMax[i] : 0);
        const lo = (i: number) => min[i] * factor + (lane ? udlMin[i] : 0);
        const points = (xs: number[], offset: number) =>
            xs.map((x, i) => ({ x, max: hi(offset + i), min: lo(offset + i) }));
        const reactionDiagrams = system.supports.map((_, s) =>
            positions.map((x, p) => ({
                x, max: reactionMax[p][s] * factor + (lane ? udlMax[system.reactionOffset + s] : 0),
                min: reactionMin[p][s] * factor + (lane ? udlMin[system.reactionOffset + s] : 0),
            })));
        const reactions = reactionDiagrams.map((diagram, s) => {
            let governing = diagram[0];
            let minimum = diagram[0].min;
            for (const point of diagram) {
                if (point.max > governing.max) governing = point;
                minimum = Math.min(minimum, point.min);
            }
            return { x: system.supports[s], max: governing.max, min: minimum, govPos: governing.x };
        });
        return {
            shear: points(system.xShear, 0), moment: points(system.xNodes, system.momentOffset),
            deflection: points(system.xNodes, system.deflectionOffset),
            reactionDiagrams, reactions, dlaUsed: lane ? 0 : dla,
        };
    };
    const cases: Partial<Record<LoadCase, CaseResults>> = {};
    if (config.loadCase !== 'lane') cases.truck = makeCase(1 + dla, false);
    if (config.loadCase !== 'truck') cases.lane = makeCase(0.8, true);
    if (cases.truck && cases.lane) {
        const truck = cases.truck;
        const lane = cases.lane;
        const combine = (t: EnvelopePoint[], l: EnvelopePoint[]) =>
            t.map((p, i) => ({ x: p.x, max: Math.max(p.max, l[i].max), min: Math.min(p.min, l[i].min) }));
        cases.envelope = {
            shear: combine(truck.shear, lane.shear),
            moment: combine(truck.moment, lane.moment),
            deflection: combine(truck.deflection, lane.deflection),
            reactionDiagrams: truck.reactionDiagrams.map((diagram, s) => combine(diagram, lane.reactionDiagrams[s])),
            reactions: truck.reactions.map((p, s) => ({
                x: p.x, max: Math.max(p.max, lane.reactions[s].max), min: Math.min(p.min, lane.reactions[s].min),
                govPos: lane.reactions[s].max > p.max ? lane.reactions[s].govPos : p.govPos,
            })),
            dlaUsed: dla,
        };
    }
    const selected = cases[config.loadCase];
    if (!selected) throw new Error('The selected load case was not calculated.');
    onProgress?.({ fraction: 1, message: 'Analysis complete.' });
    return {
        ...selected, cases, loadCase: config.loadCase, spans: spans.map(s => ({ ...s })),
        axles: axles.map(a => ({ ...a })), config: { ...config },
        xNodes: system.xNodes, supportPositions: system.supports, truckPositions: positions,
        incrementUsed: stepInfo.effective, baseIncrement: config.truckIncrement,
        incrementReason: stepInfo.reason, elapsedMs: performance.now() - started,
        stats: { factorizations: 1, truckSolves, udlSolves },
    };
}
