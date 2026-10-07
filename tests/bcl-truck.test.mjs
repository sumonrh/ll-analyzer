import test from 'node:test';
import assert from 'node:assert/strict';
import { analyzeBeam, BCL_AXLES, DEFAULT_AXLES, DEFAULT_CONFIG, buildBclSpacings,
    buildTruckPositions, computeEffectiveIncrement } from '../src/beam-engine.ts';

const spans = [12, 18, 12].map((length, i) => ({ id: `s${i}`, length }));
const close = (actual, expected) => assert.ok(Math.abs(actual - expected) <= 1e-7 * (1 + Math.abs(expected)),
    `${actual} != ${expected}`);
const config = { ...DEFAULT_CONFIG, nElemsPerSpan: 3, truckIncrement: 1, loadCase: 'envelope' };
const flatten = result => [...result.shear, ...result.moment, ...result.deflection, ...result.reactions];

test('BCL subdivisions include exact endpoints and reject unsupported values', () => {
    for (const [step, count] of [[0.5, 24], [1, 13], [2, 7]]) {
        const gaps = buildBclSpacings(step);
        assert.equal(gaps.length, count);
        assert.equal(gaps[0], 6.6);
        assert.equal(gaps.at(-1), 18);
        gaps.slice(0, -1).forEach((gap, i) => close(gap, 6.6 + step * i));
    }
    for (const value of [0, -1, 0.001, 1.5, NaN, Infinity])
        assert.throws(() => buildBclSpacings(value), /0.5, 1 or 2/);
});

test('every BCL subdivision envelopes separate fixed trucks without repeating FEM/UDL work', () => {
    for (const step of [0.5, 1, 2]) {
        const progress = [];
        const input = { spans, axles: DEFAULT_AXLES, config: { ...config, truckModel: 'BCL-625', bclSubdivision: step } };
        const snapshot = structuredClone(input);
        const actual = analyzeBeam(input, p => progress.push(p.fraction));
        assert.deepEqual(input, snapshot);
        assert.deepEqual(actual.axles, BCL_AXLES);
        assert.equal(actual.stats.factorizations, 1);
        assert.equal(actual.stats.influenceSolves, 36);
        assert.equal(actual.config.bclSubdivision, step);
        progress.forEach((p, i) => assert.ok(i === 0 || p >= progress[i - 1]));
        assert.equal(progress.at(-1), 1);
        const fixed = actual.bclSpacings.map(spacing => analyzeBeam({
            spans, axles: BCL_AXLES.map((axle, i) => i === 2 ? { ...axle, spacing } : axle),
            config: { ...config, truckModel: 'Custom' },
        }));
        assert.equal(actual.stats.udlIntegrations, fixed[0].stats.udlIntegrations);
        assert.deepEqual(actual.udlTracer, fixed[0].udlTracer);
        for (const name of ['truck', 'lane', 'envelope']) {
            flatten(actual.cases[name]).forEach((point, i) => {
                close(point.max, Math.max(...fixed.map(result => flatten(result.cases[name])[i].max)));
                close(point.min, Math.min(...fixed.map(result => flatten(result.cases[name])[i].min)));
            });
            actual.cases[name].reactions.forEach((reaction, s) => {
                assert.ok(actual.bclSpacings.includes(reaction.govSpacing));
                const governing = fixed[actual.bclSpacings.indexOf(reaction.govSpacing)].cases[name].reactions[s];
                close(reaction.max, governing.max);
                close(reaction.govPos, governing.govPos);
            });
        }
        // The common history grid spans the longest truck and includes alignments for all gaps.
        const trucks = actual.bclSpacings.map(spacing => BCL_AXLES.map((a, i) => i === 2 ? { ...a, spacing } : a));
        const longest = trucks.at(-1);
        const increment = computeEffectiveIncrement(spans, longest, config.truckIncrement, config.nElemsPerSpan).effective;
        assert.deepEqual(actual.truckPositions, buildTruckPositions(actual.supportPositions, longest, increment, trucks));
        for (const result of fixed) {
            for (const name of ['truck', 'lane']) {
                result.cases[name].reactionDiagrams.forEach((diagram, s) => {
                    diagram.forEach(point => {
                        const index = actual.truckPositions.findIndex(x => Math.abs(x - point.x) < 1e-9);
                        if (index < 0) return;
                        const envelope = actual.cases[name].reactionDiagrams[s][index];
                        assert.ok(envelope.max >= point.max - 1e-7);
                        assert.ok(envelope.min <= point.min + 1e-7);
                    });
                });
            }
        }
    }
});

test('preset normalization, Custom compatibility and manual/zero DLA lane behavior', () => {
    const custom = [{ id: 'single', load: 90, spacing: 0 }];
    const legacy = analyzeBeam({ spans, axles: custom, config });
    const explicit = analyzeBeam({ spans, axles: custom, config: { ...config, truckModel: 'Custom' } });
    assert.deepEqual(explicit.moment, legacy.moment);
    const preset = analyzeBeam({ spans, axles: custom, config: { ...config, truckModel: 'CL-625' } });
    assert.deepEqual(preset.moment, analyzeBeam({ spans, axles: DEFAULT_AXLES, config }).moment);
    assert.deepEqual(preset.bclSpacings, []);
    for (const settings of [{ dlaOverride: 0.25, dlaMultiplier: 0.75, laneUdl: 0 },
        { dlaOverride: null, dlaMultiplier: 0, laneUdl: 9 }]) {
        const result = analyzeBeam({ spans, axles: custom,
            config: { ...config, ...settings, truckModel: 'BCL-625', bclSubdivision: 2 } });
        const fixed = result.bclSpacings.map(spacing => analyzeBeam({ spans,
            axles: BCL_AXLES.map((a, i) => i === 2 ? { ...a, spacing } : a),
            config: { ...config, ...settings } }));
        for (const name of ['truck', 'lane']) flatten(result.cases[name]).forEach((point, i) => {
            close(point.max, Math.max(...fixed.map(r => flatten(r.cases[name])[i].max)));
            close(point.min, Math.min(...fixed.map(r => flatten(r.cases[name])[i].min)));
        });
    }
    assert.throws(() => analyzeBeam({ spans, axles: custom, config: { ...config, truckModel: 'BCL-625', bclSubdivision: 0.01 } }), /subdivision/);
});
