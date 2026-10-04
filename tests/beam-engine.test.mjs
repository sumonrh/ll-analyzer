import test from 'node:test';
import assert from 'node:assert/strict';
import { analyzeBeam, computeAutoDlaInfo, computeEffectiveIncrement, buildTruckPositions,
    truckGroupDla, resolveDla,
    DEFAULT_CONFIG, DEFAULT_SPANS, DEFAULT_AXLES, MAX_SWEEP_STEPS } from '../src/beam-engine.ts';
import { referenceAnalysis, referencePartialUdl } from './vba-reference.mjs';

const spansOf = lengths => lengths.map((length, i) => ({ id: `s${i}`, length }));
const axlesOf = entries => entries.map(([load, spacing], i) => ({ id: `a${i}`, load, spacing }));
const run = (spans, axles, config = {}) =>
    analyzeBeam({ spans, axles, config: { ...DEFAULT_CONFIG, ...config } });
const close = (actual, expected, tolerance = 1e-6) =>
    assert.ok(Math.abs(actual - expected) <= tolerance * Math.max(1, Math.abs(expected)),
        `${actual} != ${expected}`);
const flatten = data => [...data.shear, ...data.moment, ...data.deflection, ...data.reactions];

test('VBA parity: continuous influence optima bound the sampled VBA-equation reference', () => {
    for (const lengths of [[8], [4, 7, 3], [9, 3, 6, 4]]) {
        const spans = spansOf(lengths);
        const axles = axlesOf([[80, 1.1], [120, 2.3], [70, 0]]);
        const config = { ...DEFAULT_CONFIG, nElemsPerSpan: 6, truckIncrement: 0.2, dlaOverride: 0.3, dlaMultiplier: 1, loadCase: 'envelope' };
        const actual = run(spans, axles, config);
        const expected = referenceAnalysis(spans, axles, { ...config, step: actual.incrementUsed });
        // Influence-zone UDL (incl. partial elements) + continuous truck search must
        // retain or improve every extremum of the sampled/exhaustive reference.
        for (const name of ['truck', 'lane']) {
            flatten(actual.cases[name]).forEach((point, i) => {
                const hi = expected[name].max[i];
                const lo = expected[name].min[i];
                assert.ok(point.max >= hi - 1e-6 * Math.max(1, Math.abs(hi)), `${name}[${i}].max ${point.max} < ${hi}`);
                assert.ok(point.min <= lo + 1e-6 * Math.max(1, Math.abs(lo)), `${name}[${i}].min ${point.min} > ${lo}`);
            });
        }
    }
});

test('standard 20/25/20m CL-625 influence results bound the independent sampled reference', () => {
    const config = { ...DEFAULT_CONFIG, dlaOverride: 0.25, dlaMultiplier: 1, loadCase: 'envelope' };
    const actual = run(DEFAULT_SPANS, DEFAULT_AXLES, config);
    const expected = referenceAnalysis(DEFAULT_SPANS, DEFAULT_AXLES, { ...config, step: actual.incrementUsed });
    for (const name of ['truck', 'lane']) {
        flatten(actual.cases[name]).forEach((p, i) => {
            const hi = expected[name].max[i];
            const lo = expected[name].min[i];
            assert.ok(p.max >= hi - 2e-6 * Math.max(1, Math.abs(hi)));
            assert.ok(p.min <= lo + 2e-6 * Math.max(1, Math.abs(lo)));
        });
    }
});

test('exact-coordinate grid retains or improves every extremum of the original VBA per-pass sampling', () => {
    const spans = spansOf([4, 7, 3]);
    const axles = axlesOf([[80, 1.1], [120, 2.3], [70, 0]]);
    const config = { ...DEFAULT_CONFIG, nElemsPerSpan: 6, truckIncrement: 0.2, dlaOverride: 0.3, dlaMultiplier: 1, loadCase: 'envelope' };
    const actual = run(spans, axles, config);
    const expected = referenceAnalysis(spans, axles, { ...config, step: actual.incrementUsed, vbaSampling: true });
    for (const name of ['truck', 'lane']) {
        flatten(actual.cases[name]).forEach((p, i) => {
            const hi = expected[name].max[i];
            const lo = expected[name].min[i];
            assert.ok(p.max >= hi - 1e-6 * Math.max(1, Math.abs(hi)));
            assert.ok(p.min <= lo + 1e-6 * Math.max(1, Math.abs(lo)));
        });
    }
});

test('simply supported moving point load agrees with exact closed-form moment, deflection and reactions', () => {
    const result = run(spansOf([8]), axlesOf([[100, 0]]), { nElemsPerSpan: 8, dlaOverride: 0, dlaMultiplier: 1 });
    close(result.moment[4].max, 100 * 8 / 4);
    close(result.deflection[4].min, -100000 * 8 ** 3 / (48 * DEFAULT_CONFIG.E * DEFAULT_CONFIG.I), 1e-10);
    result.reactions.forEach(r => { close(r.max, 100); close(r.min, 0); });
    result.reactionDiagrams[0].forEach((p, i) => {
        close(p.max + result.reactionDiagrams[1][i].max, 100);
    });
    close(result.shear[0].max, 100);
    close(result.shear.at(-1).min, -100);
});

test('lane UDL agrees with exact simply supported solution; no DLA on lane truck or UDL', () => {
    const spans = spansOf([8]);
    const zero = axlesOf([[0, 0]]);
    const result = run(spans, zero, { loadCase: 'lane', nElemsPerSpan: 8, dlaOverride: 2, dlaMultiplier: 1 });
    close(result.moment[4].max, 9 * 8 ** 2 / 8);
    close(result.deflection[4].min, -5 * 9000 * 8 ** 4 / (384 * DEFAULT_CONFIG.E * DEFAULT_CONFIG.I), 1e-10);
    close(result.reactions[0].max, 9 * 8 / 2);
    assert.equal(result.dlaUsed, 0);
    const axles = axlesOf([[100, 0]]);
    const a = run(spans, axles, { loadCase: 'lane', dlaOverride: 0, dlaMultiplier: 1 });
    const b = run(spans, axles, { loadCase: 'lane', dlaOverride: 0.4, dlaMultiplier: 1 });
    assert.deepEqual(flatten(a), flatten(b));
    close(a.reactions[0].max, 80 + 36);
});

test('lane UDL intensity is configurable and defaults to 9 kN/m', () => {
    const spans = spansOf([8]);
    const zero = axlesOf([[0, 0]]);
    const nine = run(spans, zero, { loadCase: 'lane', nElemsPerSpan: 8 });
    close(nine.moment[4].max, 9 * 8 ** 2 / 8);
    close(nine.reactions[0].max, 9 * 8 / 2);
    assert.equal(nine.config.laneUdl, 9);
    const seven = run(spans, zero, { loadCase: 'lane', nElemsPerSpan: 8, laneUdl: 7 });
    close(seven.moment[4].max, 7 * 8 ** 2 / 8);
    close(seven.deflection[4].min, -5 * 7000 * 8 ** 4 / (384 * DEFAULT_CONFIG.E * DEFAULT_CONFIG.I), 1e-10);
    close(seven.reactions[0].max, 7 * 8 / 2);
    assert.equal(seven.config.laneUdl, 7);
    const eight = run(spans, zero, { loadCase: 'lane', nElemsPerSpan: 8, laneUdl: 8 });
    close(eight.moment[4].max, 8 * 8 ** 2 / 8);
});

test('DLA multiplier d scales truck DLA; auto uses selected-axle 40/30/25% groups', () => {
    assert.equal(truckGroupDla(0, false), 0);
    assert.equal(truckGroupDla(1, false), 0.4);
    assert.equal(truckGroupDla(2, false), 0.3);
    assert.equal(truckGroupDla(3, true), 0.3);
    assert.equal(truckGroupDla(3, false), 0.25);
    assert.equal(truckGroupDla(5, false), 0.25);
    assert.deepEqual(resolveDla({ ...DEFAULT_CONFIG, dlaOverride: null, dlaMultiplier: 0.75 }), { isAuto: true, base: 0, multiplier: 0.75, effective: 0 });
    assert.deepEqual(resolveDla({ ...DEFAULT_CONFIG, dlaOverride: 0.25, dlaMultiplier: 0.75 }), { isAuto: false, base: 0.25, multiplier: 0.75, effective: 0.1875 });
    const spans = spansOf([20, 25, 20]);
    const full = run(spans, DEFAULT_AXLES, { nElemsPerSpan: 8, loadCase: 'truck', dlaOverride: null, dlaMultiplier: 1 });
    const off = run(spans, DEFAULT_AXLES, { nElemsPerSpan: 8, loadCase: 'truck', dlaOverride: null, dlaMultiplier: 0 });
    const half = run(spans, DEFAULT_AXLES, { nElemsPerSpan: 8, loadCase: 'truck', dlaOverride: null, dlaMultiplier: 0.5 });
    assert.equal(full.dlaAuto, true);
    assert.equal(full.dlaMultiplier, 1);
    const maxFull = Math.max(...full.moment.map(p => p.max));
    const maxOff = Math.max(...off.moment.map(p => p.max));
    const maxHalf = Math.max(...half.moment.map(p => p.max));
    assert.ok(maxFull > maxOff, 'd=1 must exceed d=0 with auto DLA');
    assert.ok(maxHalf > maxOff && maxHalf < maxFull, 'd=0.5 must lie between off and full');
    // Override scales by d as well: 0.25 * 0.75 = 0.1875 effective.
    const over = run(spans, DEFAULT_AXLES, { nElemsPerSpan: 8, loadCase: 'truck', dlaOverride: 0.25, dlaMultiplier: 0.75 });
    assert.equal(over.dlaAuto, false);
    close(over.dlaUsed, 0.1875);
    const overFull = run(spans, DEFAULT_AXLES, { nElemsPerSpan: 8, loadCase: 'truck', dlaOverride: 0.25, dlaMultiplier: 1 });
    assert.ok(Math.max(...overFull.moment.map(p => p.max)) > Math.max(...over.moment.map(p => p.max)));
});

test('DLA evaluates every span, including short-span, tandem and coincident-axle cases', () => {
    assert.equal(computeAutoDlaInfo(spansOf([20, 0.5]), DEFAULT_AXLES).dla, 0.4);
    assert.equal(computeAutoDlaInfo(spansOf([20, 0.5]), DEFAULT_AXLES).governingSpan, 0.5);
    assert.equal(computeAutoDlaInfo(spansOf([20, 2]), DEFAULT_AXLES).dla, 0.3);
    assert.equal(computeAutoDlaInfo(DEFAULT_SPANS, DEFAULT_AXLES).dla, 0.25);
    assert.equal(computeAutoDlaInfo(spansOf([1]), axlesOf([[50, 0], [100, 0], [100, 0]])).dla, 0.25);
});

test('envelope governs V/M/D; support summaries are continuous optima bounding the sampled diagrams', () => {
    const result = run(spansOf([4, 7, 3]), DEFAULT_AXLES, { nElemsPerSpan: 8, loadCase: 'envelope' });
    const truck = flatten(result.cases.truck);
    const lane = flatten(result.cases.lane);
    flatten(result).forEach((p, i) => {
        if (i < result.shear.length + result.moment.length + result.deflection.length) {
            assert.equal(p.max, Math.max(truck[i].max, lane[i].max));
            assert.equal(p.min, Math.min(truck[i].min, lane[i].min));
        }
    });
    for (const data of Object.values(result.cases)) {
        data.reactionDiagrams.forEach((diagram, s) => {
            const summary = data.reactions[s];
            assert.ok(Number.isFinite(summary.govPos));
            assert.ok(summary.max >= Math.max(...diagram.map(p => p.max)) - 1e-9 * Math.max(1, Math.abs(summary.max)));
            assert.ok(summary.min <= Math.min(...diagram.map(p => p.min)) + 1e-9 * Math.max(1, Math.abs(summary.min)));
            assert.deepEqual(diagram.map(p => p.x), result.truckPositions);
        });
    }
    assert.equal(result.stats.factorizations, 1);
    assert.equal(result.stats.udlIntegrations, 24 * flatten(result).length);
    assert.equal(result.stats.influenceSolves, 4 * 24);
    assert.ok(result.stats.truckSolves > 0);
});

test('sweep cap and minimum floor cannot be defeated by tiny base steps; tail and exact alignments retained', () => {
    const spans = spansOf([1000, 1000]);
    const info = computeEffectiveIncrement(spans, DEFAULT_AXLES, 0.001, 40);
    const truckLength = DEFAULT_AXLES.slice(0, -1).reduce((s, a) => s + a.spacing, 0);
    assert.ok((2000 + 2 * truckLength) / info.effective <= MAX_SWEEP_STEPS);
    assert.ok(info.wasAdjusted);
    assert.ok(!info.wasReduced);
    const points = buildTruckPositions([0, 1000, 2000], DEFAULT_AXLES, info.effective);
    assert.equal(points[0], -truckLength);
    assert.equal(points.at(-1), 2000 + truckLength);
    assert.ok(points.includes(1000 + 3.6));
    assert.ok(points.includes(1000 - 3.6));
    assert.equal(computeEffectiveIncrement(spansOf([1]), axlesOf([[50, 0]]), 0.001, 40).effective, 0.02);
    assert.equal(computeEffectiveIncrement(spansOf([8]), axlesOf([[50, 0]]), 0.023, 40).effective, 0.023);
});

test('output shapes are stepped shear, nodal moment/deflection and reaction histories at exact coordinates', () => {
    const result = run(spansOf([0.7, 4.3, 1.1]), DEFAULT_AXLES, { nElemsPerSpan: 10, loadCase: 'envelope' });
    assert.equal(result.shear.length, 60);
    assert.equal(result.moment.length, 31);
    assert.equal(result.deflection.length, 31);
    assert.equal(result.reactions.length, 4);
    assert.equal(result.shear[1].x, result.shear[2].x);
    for (const data of Object.values(result.cases)) {
        for (const p of flatten(data)) {
            assert.ok(Number.isFinite(p.max) && Number.isFinite(p.min) && p.max >= p.min);
        }
        for (let s = 0; s < 4; s++) {
            close(data.deflection[s * 10].max, 0, 1e-12);
            close(data.deflection[s * 10].min, 0, 1e-12);
            assert.deepEqual(data.reactionDiagrams[s].map(p => p.x), result.truckPositions);
        }
    }
});

test('inputs are never mutated and the unused trailing axle spacing has no effect', () => {
    const spans = structuredClone(DEFAULT_SPANS);
    const axles = structuredClone(DEFAULT_AXLES);
    const before = JSON.stringify({ spans, axles });
    const a = run(spans, axles, { nElemsPerSpan: 4 });
    assert.equal(JSON.stringify({ spans, axles }), before);
    axles.at(-1).spacing = 999;
    const b = run(spans, axles, { nElemsPerSpan: 4 });
    assert.deepEqual(flatten(a), flatten(b));
    assert.deepEqual(a.truckPositions, b.truckPositions);
});

test('EI scaling preserves force results and scales deflections', () => {
    const a = run(spansOf([5, 9]), DEFAULT_AXLES, { nElemsPerSpan: 6 });
    const b = run(spansOf([5, 9]), DEFAULT_AXLES, { nElemsPerSpan: 6, E: DEFAULT_CONFIG.E * 2 });
    a.moment.forEach((p, i) => { close(p.max, b.moment[i].max); close(p.min, b.moment[i].min); });
    a.deflection.forEach((p, i) => { close(p.max / 2, b.deflection[i].max, 1e-12); close(p.min / 2, b.deflection[i].min, 1e-12); });
});

test('reversing the axle input preserves static (d=0) envelopes on unequal spans', () => {
    const reversed = DEFAULT_AXLES.map((_, i) => ({
        ...DEFAULT_AXLES[DEFAULT_AXLES.length - 1 - i],
        spacing: i < DEFAULT_AXLES.length - 1 ? DEFAULT_AXLES[DEFAULT_AXLES.length - 2 - i].spacing : 0,
    }));
    const opts = { nElemsPerSpan: 6, loadCase: 'envelope', dlaMultiplier: 0 };
    const a = run(spansOf([3, 7, 5]), DEFAULT_AXLES, opts);
    const b = run(spansOf([3, 7, 5]), reversed, opts);
    flatten(a).forEach((p, i) => {
        close(p.max, flatten(b)[i].max);
        close(p.min, flatten(b)[i].min);
    });
});

test('coincident axles and loads exactly at supports retain equilibrium', () => {
    const result = run(spansOf([8]), axlesOf([[50, 0], [75, 0]]), { nElemsPerSpan: 8, dlaOverride: 0, dlaMultiplier: 1 });
    close(result.moment[4].max, 125 * 8 / 4);
    result.reactionDiagrams[0].forEach((point, p) => {
        close(point.max + result.reactionDiagrams[1][p].max, 125);
    });
    assert.equal(result.reactions[0].govPos, 0);
    assert.equal(result.reactions[1].govPos, 8);
});

test('invalid input produces explicit errors before allocation or solving', () => {
    for (const config of [{ nElemsPerSpan: 2.5 }, { nElemsPerSpan: 1 }, { truckIncrement: 0 },
        { E: NaN }, { I: -1 }, { dlaOverride: -0.1 }, { dlaMultiplier: -0.1 }, { dlaMultiplier: 1.5 }, { laneUdl: -1 }, { laneUdl: NaN }, { loadCase: 'invalid' }])
        assert.throws(() => run(DEFAULT_SPANS, DEFAULT_AXLES, config));
    assert.throws(() => run([], DEFAULT_AXLES), /span/);
    assert.throws(() => run(spansOf([0]), DEFAULT_AXLES), /Span 1/);
    assert.throws(() => run(DEFAULT_SPANS, []), /axles/);
    assert.throws(() => run(DEFAULT_SPANS, axlesOf([[-1, 0]])), /Axle 1/);
    assert.throws(() => run(DEFAULT_SPANS, axlesOf([[50, -1], [50, 0]])), /spacing/);
    assert.throws(() => run(spansOf(Array(13).fill(10)), DEFAULT_AXLES, { loadCase: 'lane' }), /12 spans/);
    assert.throws(() => run(DEFAULT_SPANS, DEFAULT_AXLES, { nElemsPerSpan: 10001 }), /1000/);
    assert.throws(() => run(spansOf(Array(400).fill(10)), axlesOf([[100, 0]]),
        { nElemsPerSpan: 2 }), /2,000,000/);
});

test('12-span lane patterning uses influence zones and finishes within 15 seconds', t => {
    const result = run(spansOf(Array(12).fill(10)), DEFAULT_AXLES, { loadCase: 'envelope' });
    assert.equal(result.stats.udlIntegrations, 480 * flatten(result).length);
    assert.equal(result.stats.factorizations, 1);
    assert.ok(result.stats.truckSolves > 0);
    assert.ok(result.elapsedMs < 15000, `Analysis took ${result.elapsedMs.toFixed(0)}ms`);
    t.diagnostic(`12-span / 40-element envelope: ${result.elapsedMs.toFixed(1)}ms, ${result.stats.truckSolves} truck optimisations`);
});

test('continuous support uplift matches the exact two-span solution, not sampled minima', () => {
    const spans = spansOf([8, 8]);
    const result = run(spans, axlesOf([[100, 0]]), {
        nElemsPerSpan: 4, dlaOverride: 0, loadCase: 'envelope',
    });
    const truckMin = -100 / (6 * Math.sqrt(3));
    const laneMin = 0.8 * truckMin - 9 * 8 / 16;
    for (const s of [0, 2]) {
        close(result.cases.truck.reactions[s].min, truckMin, 1e-10);
        close(result.cases.lane.reactions[s].min, laneMin, 1e-10);
        close(result.reactions[s].min, laneMin, 1e-10);
        const sampled = Math.min(...result.cases.truck.reactionDiagrams[s].map(p => p.min));
        assert.ok(sampled - result.cases.truck.reactions[s].min > 1e-5,
            'This benchmark must distinguish the sampled and continuous minimum');
    }
});

test('blank lane UDL defaults to 9; zero preserves the 80% lane truck and disables only UDL', () => {
    const spans = spansOf([4, 7, 3]);
    const options = { nElemsPerSpan: 6, loadCase: 'lane' };
    const blank = run(spans, DEFAULT_AXLES, { ...options, laneUdl: null });
    const omitted = run(spans, DEFAULT_AXLES, { ...options, laneUdl: undefined });
    const nine = run(spans, DEFAULT_AXLES, { ...options, laneUdl: 9 });
    assert.equal(blank.config.laneUdl, 9);
    assert.deepEqual(flatten(blank), flatten(nine));
    assert.deepEqual(flatten(omitted), flatten(nine));
    const zero = run(spans, DEFAULT_AXLES, { ...options, laneUdl: 0, dlaOverride: 0.9 });
    const staticTruck = run(spans, DEFAULT_AXLES, { nElemsPerSpan: 6, dlaOverride: 0 });
    flatten(zero).forEach((p, i) => {
        close(p.max, 0.8 * flatten(staticTruck)[i].max, 1e-10);
        close(p.min, 0.8 * flatten(staticTruck)[i].min, 1e-10);
    });
    assert.equal(zero.config.laneUdl, 0);
    assert.equal(zero.stats.udlIntegrations, 0);
    assert.equal(zero.udlTracer.intensity, 0);
    assert.equal(zero.udlTracer.max, 0);
    assert.equal(zero.udlTracer.min, 0);
    close(zero.udlTracer.reconstructedMax, 0);
    close(zero.udlTracer.reconstructedMin, 0);
    assert.deepEqual(zero.udlTracer.intervals, []);
    assert.deepEqual(zero.udlTracer.influenceLine, nine.udlTracer.influenceLine);
    assert.ok(zero.udlTracer.influenceLine.some(p => Math.abs(p.ordinate) > 0.1));
});

test('automatic UDL tracer selects the nearest midpoint node and verifies exact partial-element loading', () => {
    for (const lengths of [[8], [2, 15, 7], [4, 9, 1]]) {
        const spans = spansOf(lengths);
        const config = { ...DEFAULT_CONFIG, nElemsPerSpan: 4, loadCase: 'lane', laneUdl: 7 };
        const result = run(spans, axlesOf([[0, 0]]), config);
        const trace = result.udlTracer;
        const midpoint = lengths.reduce((a, b) => a + b) / 2;
        const nearest = result.xNodes.reduce((best, x, i) =>
            Math.abs(x - midpoint) < Math.abs(result.xNodes[best] - midpoint) ? i : best, 0);
        assert.equal(trace.nodeIndex, nearest + 1);
        assert.equal(trace.x, result.xNodes[nearest]);
        close(trace.max, result.moment[nearest].max, 1e-10);
        close(trace.min, result.moment[nearest].min, 1e-10);
        for (const envelope of ['max', 'min']) {
            const intervals = trace.intervals.filter(i => i.envelope === envelope);
            intervals.forEach(i => {
                assert.ok(i.xiStart >= 0 && i.xiEnd <= 1 && i.xiEnd > i.xiStart);
                assert.ok(envelope === 'max' ? i.contribution > 0 : i.contribution < 0);
            });
            close(intervals.reduce((sum, i) => sum + i.contribution, 0), trace[envelope], 1e-10);
            const independent = referencePartialUdl(spans, config, intervals, trace.intensity);
            close(independent.moment[nearest], trace[envelope], 1e-9);
        }
        assert.ok(Math.abs(trace.reconstructedMax - trace.max) <= trace.tolerance);
        assert.ok(Math.abs(trace.reconstructedMin - trace.min) <= trace.tolerance);
        if (lengths.length > 1) {
            assert.ok(trace.intervals.some(i => i.xiStart > 1e-8 || i.xiEnd < 1 - 1e-8),
                'The benchmark must exercise a partial element, not just whole-span patterning');
        }
    }
    const truck = run(spansOf([8]), axlesOf([[100, 0]]), { nElemsPerSpan: 4 });
    assert.equal(truck.udlTracer, undefined);
});

test('unit influence solves stay in the intended element for very short elements', () => {
    const result = run(spansOf([8e-6]), axlesOf([[0, 0]]), {
        loadCase: 'lane', nElemsPerSpan: 16,
    });
    close(result.udlTracer.max, 9 * (8e-6) ** 2 / 8, 1e-20);
    close(result.udlTracer.reconstructedMax, result.udlTracer.max, 1e-20);
});

test('progress remains monotonic through both cases and UDL verification', () => {
    const progress = [];
    analyzeBeam({
        spans: spansOf([4, 7, 3]), axles: DEFAULT_AXLES,
        config: { ...DEFAULT_CONFIG, nElemsPerSpan: 4, loadCase: 'envelope' },
    }, p => progress.push(p.fraction));
    progress.forEach((f, i) => assert.ok(f >= 0 && f <= 1 && (i === 0 || f >= progress[i - 1])));
    assert.equal(progress.at(-1), 1);
});
