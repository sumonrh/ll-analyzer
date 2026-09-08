import { execFileSync } from 'node:child_process';
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = dirname(fileURLToPath(import.meta.url));
const appDir = resolve(root, 'll-analyzer');
const distDir = resolve(appDir, 'dist');
const assetsDir = resolve(distDir, 'assets');
const rootIndexPath = resolve(root, 'index.html');
const indexPath = resolve(appDir, 'index.html');
const distIndexPath = resolve(distDir, 'index.html');
const templatePath = resolve(appDir, 'index.template.html');
const npmCommand = process.platform === 'win32' ? 'npm.cmd' : 'npm';

const existingIndex = readFileSync(indexPath, 'utf8');
const buildTemplate = readFileSync(templatePath, 'utf8');

try {
  writeFileSync(indexPath, buildTemplate, 'utf8');

  execFileSync(npmCommand, ['run', 'build:app', '--workspace=ll-analyzer'], {
    cwd: root,
    shell: process.platform === 'win32',
    stdio: 'inherit',
  });

  const builtHtml = readFileSync(resolve(distDir, 'index.html'), 'utf8');
  const jsFile = builtHtml.match(/<script[^>]+src="\.\/assets\/([^"']+\.js)"/)?.[1];
  const cssFile = builtHtml.match(/<link[^>]+href="\.\/assets\/([^"']+\.css)"/)?.[1];

  if (!jsFile || !cssFile) {
    throw new Error('Could not find the built JavaScript and CSS assets.');
  }

  const jsContent = readFileSync(resolve(assetsDir, jsFile), 'utf8')
    .replace(/<\/script/gi, '<\\/script');
  const cssContent = readFileSync(resolve(assetsDir, cssFile), 'utf8')
    .replace(/<\/style/gi, '<\\/style');

  const standaloneHtml = `<!doctype html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>LL Analyzer</title>
    <style>
${cssContent}
    </style>
  </head>
  <body>
    <div id="root"></div>
    <script type="module">
${jsContent}
    </script>
  </body>
</html>
`;

  writeFileSync(rootIndexPath, standaloneHtml, 'utf8');
  writeFileSync(indexPath, standaloneHtml, 'utf8');
  writeFileSync(distIndexPath, standaloneHtml, 'utf8');
  console.log(`Successfully generated standalone index.html at ${rootIndexPath} and ${indexPath} (${(standaloneHtml.length / 1024).toFixed(0)} KB)`);
} catch (error) {
  writeFileSync(indexPath, existingIndex, 'utf8');
  throw error;
}