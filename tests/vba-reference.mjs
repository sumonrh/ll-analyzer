// Independent dense FEM implementation; the legacy sweep enumerates whole-span UDL patterns.
function referenceSystem(spans, config) {
    const n = spans.length * config.nElemsPerSpan;
    const nd = 2 * (n + 1);
    const x = [0];
    const lengths = [];
    const supports = [0];
    for (const span of spans) {
        const start = supports.at(-1);
        const le = span.length / config.nElemsPerSpan;
        for (let i = 1; i <= config.nElemsPerSpan; i++) {
            x.push(start + i * le);
            lengths.push(le);
        }
        supports.push(start + span.length);
    }
    const local = lengths.map(le => {
        const c = config.E * config.I / le ** 3;
        return [
            [12, 6 * le, -12, 6 * le],
            [6 * le, 4 * le * le, -6 * le, 2 * le * le],
            [-12, -6 * le, 12, -6 * le],
            [6 * le, 2 * le * le, -6 * le, 4 * le * le],
        ].map(row => row.map(v => v * c));
    });
    const K = Array.from({ length: nd }, () => Array(nd).fill(0));
    local.forEach((matrix, e) => {
        for (let r = 0; r < 4; r++)
            for (let c = 0; c < 4; c++) K[e * 2 + r][e * 2 + c] += matrix[r][c];
    });
    const constrained = supports.map((_, s) => s * config.nElemsPerSpan * 2);
    const free = Array.from({ length: nd }, (_, i) => i).filter(i => !constrained.includes(i));
    const LU = free.map(i => free.map(j => K[i][j]));
    for (let k = 0; k < free.length - 1; k++) {
        for (let i = k + 1; i < free.length; i++) {
            LU[i][k] /= LU[k][k];
            for (let j = k + 1; j < free.length; j++) LU[i][j] -= LU[i][k] * LU[k][j];
        }
    }
    const solve = (F, elemLoads) => {
        const y = free.map(i => F[i]);
        for (let i = 1; i < free.length; i++)
            for (let j = 0; j < i; j++) y[i] -= LU[i][j] * y[j];
        for (let i = free.length - 1; i >= 0; i--) {
            for (let j = i + 1; j < free.length; j++) y[i] -= LU[i][j] * y[j];
            y[i] /= LU[i][i];
        }
        const U = Array(nd).fill(0);
        free.forEach((dof, i) => { U[dof] = y[i]; });
        const forces = local.map((matrix, e) => matrix.map((row, r) =>
            row.reduce((sum, v, c) => sum + v * U[e * 2 + c], 0) - elemLoads[e][r]));
        const shear = forces.flatMap(f => [f[0] / 1000, -f[2] / 1000]);
        const moment = x.map((_, node) => {
            const value = node === 0 ? forces[0][1] : node === n ? -forces[n - 1][3]
                : (-forces[node - 1][3] + forces[node][1]) / 2;
            return -value / 1000;
        });
        const deflection = x.map((_, node) => U[node * 2]);
        const reactions = constrained.map(dof =>
            (K[dof].reduce((sum, k, j) => sum + k * U[j], 0) - F[dof]) / 1000);
        return [...shear, ...moment, ...deflection, ...reactions];
    };
    const loadArrays = () => ({
        F: Array(nd).fill(0), elemLoads: Array.from({ length: n }, () => [0, 0, 0, 0]),
    });
    return { n, x, lengths, supports, solve, loadArrays };
}

export function referencePartialUdl(spans, config, intervals, intensity) {
    const { n, x, lengths, solve, loadArrays } = referenceSystem(spans, config);
    const { F, elemLoads } = loadArrays();
    // Two-point Gauss quadrature integrates the cubic Hermite shapes exactly.
    for (const interval of intervals) {
        const e = interval.element - 1;
        const le = lengths[e];
        const a = (interval.start - x[e]) / le;
        const b = (interval.end - x[e]) / le;
        const half = (b - a) / 2;
        for (const sign of [-1, 1]) {
            const t = (a + b) / 2 + sign * half / Math.sqrt(3);
            const shapes = [1 - 3 * t ** 2 + 2 * t ** 3, le * (t - 2 * t ** 2 + t ** 3),
                3 * t ** 2 - 2 * t ** 3, le * (-(t ** 2) + t ** 3)];
            for (let r = 0; r < 4; r++) {
                const value = -intensity * 1000 * le * half * shapes[r];
                F[e * 2 + r] += value;
                elemLoads[e][r] += value;
            }
        }
    }
    const packed = solve(F, elemLoads);
    return { moment: packed.slice(2 * n, 3 * n + 1) };
}

export function referenceAnalysis(spans, axles, config) {
    const { n, x, lengths, supports, solve, loadArrays } = referenceSystem(spans, config);
    const truckResponse = (lead, direction) => {
        const { F, elemLoads } = loadArrays();
        let pos = lead;
        for (let a = 0; a < direction.length; a++) {
            if (pos >= -1e-4 && pos <= supports.at(-1) + 1e-4) {
                const bounded = Math.min(supports.at(-1), Math.max(0, pos));
                let e = lengths.findIndex((_, i) => bounded <= x[i + 1] + 1e-6);
                if (e < 0) e = n - 1;
                const le = lengths[e];
                const t = Math.min(1, Math.max(0, (bounded - x[e]) / le));
                const shapes = [1 - 3 * t ** 2 + 2 * t ** 3, le * (t - 2 * t ** 2 + t ** 3),
                    3 * t ** 2 - 2 * t ** 3, le * (-(t ** 2) + t ** 3)];
                for (let r = 0; r < 4; r++) {
                    const value = -direction[a].load * 1000 * shapes[r];
                    F[e * 2 + r] += value;
                    elemLoads[e][r] += value;
                }
            }
            if (a < direction.length - 1) pos -= direction[a].spacing;
        }
        return solve(F, elemLoads);
    };
    const size = n * 2 + 2 * (n + 1) + supports.length;
    const udlMax = Array(size).fill(0);
    const udlMin = Array(size).fill(0);
    if (config.loadCase !== 'truck') {
        for (let pattern = 0; pattern < 2 ** spans.length; pattern++) {
            const { F, elemLoads } = loadArrays();
            for (let e = 0; e < n; e++) {
                const s = Math.floor(e / config.nElemsPerSpan);
                if ((pattern & (2 ** s)) === 0) continue;
                const le = lengths[e];
                const w = config.laneUdl ?? 9;
                const values = [-w * le * 1000 / 2, -w * le ** 2 * 1000 / 12,
                    -w * le * 1000 / 2, w * le ** 2 * 1000 / 12];
                for (let r = 0; r < 4; r++) {
                    F[e * 2 + r] += values[r];
                    elemLoads[e][r] += values[r];
                }
            }
            solve(F, elemLoads).forEach((v, i) => {
                udlMax[i] = Math.max(udlMax[i], v);
                udlMin[i] = Math.min(udlMin[i], v);
            });
        }
    }
    const reversed = axles.map((_, i) => ({
        ...axles[axles.length - 1 - i],
        spacing: i < axles.length - 1 ? axles[axles.length - 2 - i].spacing : 0,
    }));
    const truckLen = axles.slice(0, -1).reduce((sum, a) => sum + a.spacing, 0);
    const start = -truckLen;
    const end = supports.at(-1) + truckLen;
    const positions = [];
    const passPositions = [];
    const count = Math.max(2, Math.floor((end - start) / config.step + 1e-6) + 1);
    for (let i = 0; i < count; i++) positions.push(Math.min(end, start + i * config.step));
    positions[positions.length - 1] = end;
    for (const direction of [axles, reversed]) {
        const samples = positions.slice(0, count);
        for (const support of supports) {
            let dist = 0;
            for (let i = 0; i < direction.length; i++) {
                for (const sign of [1, -1]) {
                    const lead = support + sign * dist;
                    if (lead >= start && lead <= end) {
                        positions.push(lead);
                        samples.push(lead);
                    }
                }
                passPositions.push(samples);
                if (i < direction.length - 1) dist += direction[i].spacing;
            }
        }
    }
    const sorted = [...new Set(positions)].sort((a, b) => a - b)
        .filter((v, i, all) => i === 0 || v - all[i - 1] > 1e-9);
    const runCase = (factor, lane) => {
        const max = Array(size).fill(-Infinity);
        const min = Array(size).fill(Infinity);
        const histories = sorted.map(() => ({ max: Array(supports.length).fill(-Infinity), min: Array(supports.length).fill(Infinity) }));
        for (const [directionIndex, direction] of [axles, reversed].entries()) {
            const factored = direction.map(a => ({ ...a, load: a.load * factor }));
            const samples = config.vbaSampling ? passPositions[directionIndex] : sorted;
            samples.forEach((lead, p) => {
                const values = truckResponse(lead, factored);
                values.forEach((value, i) => {
                    const hi = value + (lane ? udlMax[i] : 0);
                    const lo = value + (lane ? udlMin[i] : 0);
                    max[i] = Math.max(max[i], hi);
                    min[i] = Math.min(min[i], lo);
                    if (!config.vbaSampling && i >= size - supports.length) {
                        const s = i - (size - supports.length);
                        histories[p].max[s] = Math.max(histories[p].max[s], hi);
                        histories[p].min[s] = Math.min(histories[p].min[s], lo);
                    }
                });
            });
        }
        return { max, min, histories };
    };
    return {
        positions: sorted,
        truck: config.loadCase !== 'lane' ? runCase(1 + config.dlaOverride, false) : null,
        lane: config.loadCase !== 'truck' ? runCase(0.8, true) : null,
    };
}
