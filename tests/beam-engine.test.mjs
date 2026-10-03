import test from 'node:test';
import assert from 'node:assert/strict';
import { analyzeBeam, computeAutoDlaInfo, computeEffectiveIncrement, buildTruckPositions,
    DEFAULT_CONFIG, DEFAULT_SPANS, DEFAULT_AXLES, MAX_SWEEP_STEPS } from '../src/beam-engine.ts';
import { referenceAnalysis } from './vba-reference.mjs';

const spansOf = lengths => lengths.map((length, i) => ({ id: `s${i}`, length }));
const axlesOf = entries => entries.map(([load, spacing], i) => ({ id: `a${i}`, load, spacing }));
const run = (spans, axles, config = {}) =>
    analyzeBeam({ spans, axles, config: { ...DEFAULT_CONFIG, ...config } });
const close = (actual, expected, tolerance = 1e-6) =>
    assert.ok(Math.abs(actual - expected) <= tolerance * Math.max(1, Math.abs(expected)),
        `${actual} != ${expected}`);
const flatten = data => [...data.shear, ...data.moment, ...data.deflection, ...data.reactions];

test('VBA equations: banded solve and linear patterning match dense LU/exhaustive patterns', () => {
    for (const lengths of [[8], [4, 7, 3], [9, 3, 6, 4]]) {
        const spans = spansOf(lengths);
        const axles = axlesOf([[80, 1.1], [120, 2.3], [70, 0]]);
        const config = { ...DEFAULT_CONFIG, nElemsPerSpan: 6, truckIncrement: 0.2, dlaOverride: 0.3, loadCase: 'envelope' };
        const actual = run(spans, axles, config);
        const expected = referenceAnalysis(spans, axles, { ...config, step: actual.incrementUsed });
        assert.deepEqual(actual.truckPositions, expected.positions);
        for (const name of ['truck', 'lane']) {
            flatten(actual.cases[name]).forEach((point, i) => {
                close(point.max, expected[name].max[i]);
                close(point.min, expected[name].min[i]);
            });
            actual.cases[name].reactionDiagrams.forEach((diagram, s) => {
                diagram.forEach((point, p) => {
                    close(point.max, expected[name].histories[p].max[s]);
                    close(point.min, expected[name].histories[p].min[s]);
                });
            });
        }
    }
});

test('standard 20/25/20m CL-625 cases match the independent VBA-equation reference', () => {
    const config = { ...DEFAULT_CONFIG, dlaOverride: 0.25, loadCase: 'envelope' };
    const actual = run(DEFAULT_SPANS, DEFAULT_AXLES, config);
    const expected = referenceAnalysis(DEFAULT_SPANS, DEFAULT_AXLES, { ...config, step: actual.incrementUsed });
    for (const name of ['truck', 'lane']) {
        flatten(actual.cases[name]).forEach((p, i) => {
            close(p.max, expected[name].max[i], 2e-6);
            close(p.min, expected[name].min[i], 2e-6);
        });
    }
});

test('exact-coordinate grid retains or improves every extremum of the original VBA per-pass sampling', () => {
    const spans = spansOf([4, 7, 3]);
    const axles = axlesOf([[80, 1.1], [120, 2.3], [70, 0]]);
    const config = { ...DEFAULT_CONFIG, nElemsPerSpan: 6, truckIncrement: 0.2, dlaOverride: 0.3, loadCase: 'envelope' };
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
    const result = run(spansOf([8]), axlesOf([[100, 0]]), { nElemsPerSpan: 8, dlaOverride: 0 });
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
    const result = run(spans, zero, { loadCase: 'lane', nElemsPerSpan: 8, dlaOverride: 2 });
    close(result.moment[4].max, 9 * 8 ** 2 / 8);
    close(result.deflection[4].min, -5 * 9000 * 8 ** 4 / (384 * DEFAULT_CONFIG.E * DEFAULT_CONFIG.I), 1e-10);
    close(result.reactions[0].max, 9 * 8 / 2);
    assert.equal(result.dlaUsed, 0);
    const axles = axlesOf([[100, 0]]);
    const a = run(spans, axles, { loadCase: 'lane', dlaOverride: 0 });
    const b = run(spans, axles, { loadCase: 'lane', dlaOverride: 0.4 });
    assert.deepEqual(flatten(a), flatten(b));
    close(a.reactions[0].max, 80 + 36);
});

test('DLA evaluates every span, including short-span, tandem and coincident-axle cases', () => {
    assert.equal(computeAutoDlaInfo(spansOf([20, 0.5]), DEFAULT_AXLES).dla, 0.4);
    assert.equal(computeAutoDlaInfo(spansOf([20, 0.5]), DEFAULT_AXLES).governingSpan, 0.5);
    assert.equal(computeAutoDlaInfo(spansOf([20, 2]), DEFAULT_AXLES).dla, 0.3);
    assert.equal(computeAutoDlaInfo(DEFAULT_SPANS, DEFAULT_AXLES).dla, 0.25);
    assert.equal(computeAutoDlaInfo(spansOf([1]), axlesOf([[50, 0], [100, 0], [100, 0]])).dla, 0.25);
});

test('envelope includes separate cases and exactly governs every result and reaction ordinate', () => {
    const result = run(spansOf([4, 7, 3]), DEFAULT_AXLES, { nElemsPerSpan: 8, loadCase: 'envelope' });
    const truck = flatten(result.cases.truck);
    const lane = flatten(result.cases.lane);
    flatten(result).forEach((p, i) => {
        assert.equal(p.max, Math.max(truck[i].max, lane[i].max));
        assert.equal(p.min, Math.min(truck[i].min, lane[i].min));
    });
    for (const data of Object.values(result.cases)) {
        data.reactionDiagrams.forEach((diagram, s) => {
            const governing = diagram.find(p => p.x === data.reactions[s].govPos);
            assert.equal(governing.max, data.reactions[s].max);
            assert.equal(Math.min(...diagram.map(p => p.min)), data.reactions[s].min);
        });
    }
    assert.equal(result.stats.factorizations, 1);
    assert.equal(result.stats.udlSolves, 3);
    assert.equal(result.stats.truckSolves, result.truckPositions.length * 2);
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

test('reversing the axle input preserves both-direction envelopes on unequal spans', () => {
    const reversed = DEFAULT_AXLES.map((_, i) => ({
        ...DEFAULT_AXLES[DEFAULT_AXLES.length - 1 - i],
        spacing: i < DEFAULT_AXLES.length - 1 ? DEFAULT_AXLES[DEFAULT_AXLES.length - 2 - i].spacing : 0,
    }));
    const a = run(spansOf([3, 7, 5]), DEFAULT_AXLES, { nElemsPerSpan: 6, loadCase: 'envelope' });
    const b = run(spansOf([3, 7, 5]), reversed, { nElemsPerSpan: 6, loadCase: 'envelope' });
    flatten(a).forEach((p, i) => {
        close(p.max, flatten(b)[i].max);
        close(p.min, flatten(b)[i].min);
    });
});

test('coincident axles and loads exactly at supports retain equilibrium', () => {
    const result = run(spansOf([8]), axlesOf([[50, 0], [75, 0]]), { nElemsPerSpan: 8, dlaOverride: 0 });
    close(result.moment[4].max, 125 * 8 / 4);
    result.reactionDiagrams[0].forEach((point, p) => {
        close(point.max + result.reactionDiagrams[1][p].max, 125);
    });
    assert.equal(result.reactions[0].govPos, 0);
    assert.equal(result.reactions[1].govPos, 8);
});

test('invalid input produces explicit errors before allocation or solving', () => {
    for (const config of [{ nElemsPerSpan: 2.5 }, { nElemsPerSpan: 1 }, { truckIncrement: 0 },
        { E: NaN }, { I: -1 }, { dlaOverride: -0.1 }, { loadCase: 'invalid' }])
        assert.throws(() => run(DEFAULT_SPANS, DEFAULT_AXLES, config));
    assert.throws(() => run([], DEFAULT_AXLES), /span/);
    assert.throws(() => run(spansOf([0]), DEFAULT_AXLES), /Span 1/);
    assert.throws(() => run(DEFAULT_SPANS, []), /axles/);
    assert.throws(() => run(DEFAULT_SPANS, axlesOf([[-1, 0]])), /Axle 1/);
    assert.throws(() => run(DEFAULT_SPANS, axlesOf([[50, -1], [50, 0]])), /spacing/);
    assert.throws(() => run(spansOf(Array(13).fill(10)), DEFAULT_AXLES, { loadCase: 'lane' }), /12 spans/);
    assert.throws(() => run(DEFAULT_SPANS, DEFAULT_AXLES, { nElemsPerSpan: 10001 }), /10,000/);
    assert.throws(() => run(spansOf(Array(1000).fill(10)), axlesOf([[100, 0]]),
        { nElemsPerSpan: 2 }), /2,000,000/);
});

test('12-span lane patterning uses 12 rather than 4096 UDL solves and finishes within 5 seconds', t => {
    const result = run(spansOf(Array(12).fill(10)), DEFAULT_AXLES, { loadCase: 'envelope' });
    assert.equal(result.stats.udlSolves, 12);
    assert.equal(result.stats.factorizations, 1);
    assert.ok(result.elapsedMs < 5000, `Analysis took ${result.elapsedMs.toFixed(0)}ms`);
    t.diagnostic(`12-span / 40-element envelope: ${result.elapsedMs.toFixed(1)}ms, ${result.stats.truckSolves} truck solves`);
});
