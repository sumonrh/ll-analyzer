import test from 'node:test';
import assert from 'node:assert/strict';
import { copyFileSync, existsSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { spawnSync } from 'node:child_process';
import { tmpdir } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const readJson = name => JSON.parse(readFileSync(join(root, name), 'utf8'));

test('source, tests and build/deployment configuration live directly in the root', () => {
    for (const file of [join('src', 'main.tsx'), join('src', 'beam-engine.ts'), join('src', 'analysis.worker.ts'),
        'index.template.html', 'build-standalone.mjs', 'vite.config.ts', 'eslint.config.js',
        'tsconfig.json', 'tsconfig.app.json', 'tsconfig.node.json', 'firebase.json', '.firebaserc']) {
        assert.ok(existsSync(join(root, file)), `Missing root-level project file: ${file}`);
    }
    assert.ok(!existsSync(join(root, 'll-analyzer')), 'Obsolete nested app directory remains');
    const manifest = readJson('package.json');
    assert.equal(manifest.workspaces, undefined);
    assert.equal(manifest.type, 'module');
    assert.equal(manifest.scripts.build, 'node build-standalone.mjs');
    assert.equal(manifest.scripts.dev, 'vite');
    assert.equal(manifest.scripts.start, 'serve -s dist -l 8080');
    for (const script of Object.values(manifest.scripts))
        assert.ok(!script.includes('--workspace') && !script.includes('../build-standalone'));
    assert.match(readFileSync(join(root, 'apphosting.yaml'), 'utf8'), /^rootDirectory: \.$/m);
    assert.equal(readJson('firebase.json').hosting.public, 'dist');
});

test('root dependency lock matches the manifest and contains no workspace links', () => {
    const manifest = readJson('package.json');
    const lock = readJson('package-lock.json');
    assert.equal(lock.name, manifest.name);
    assert.equal(lock.version, manifest.version);
    assert.deepEqual(lock.packages[''].dependencies, manifest.dependencies);
    assert.deepEqual(lock.packages[''].optionalDependencies, manifest.optionalDependencies);
    assert.equal(lock.packages[''].workspaces, undefined);
    for (const [name, entry] of Object.entries(lock.packages)) {
        assert.ok(name === '' || name.startsWith('node_modules/'), `Obsolete project entry: ${name}`);
        assert.ok(!entry.link, `Obsolete workspace link: ${name}`);
    }
});

test('a failed root build restores the existing standalone app', () => {
    const fixture = mkdtempSync(join(tmpdir(), 'll-analyzer-flat-build-'));
    try {
        copyFileSync(join(root, 'build-standalone.mjs'), join(fixture, 'build-standalone.mjs'));
        const existing = '<!doctype html><title>Existing standalone app</title>';
        writeFileSync(join(fixture, 'index.html'), existing);
        writeFileSync(join(fixture, 'index.template.html'), '<!doctype html><title>Build template</title>');
        writeFileSync(join(fixture, 'package.json'), JSON.stringify({
            name: 'build-failure-fixture', private: true, type: 'module',
            scripts: { 'build:app': 'node -e "process.stderr.write(\'EXPECTED_BUILD_FAILURE\');process.exit(1)"' },
        }));
        const result = spawnSync(process.execPath, [join(fixture, 'build-standalone.mjs')], {
            cwd: fixture, encoding: 'utf8', timeout: 30000,
        });
        assert.ifError(result.error);
        assert.notEqual(result.status, 0);
        assert.match(result.stderr, /EXPECTED_BUILD_FAILURE/);
        assert.equal(readFileSync(join(fixture, 'index.html'), 'utf8'), existing);
    } finally {
        rmSync(fixture, { recursive: true, force: true });
    }
});
