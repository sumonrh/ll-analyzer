# LL Analyzer

LL Analyzer is a modern, professional web application designed for comprehensive Live Load Analysis of continuous beams and out-of-the-box NU Girder designs. Built with React, Vite, and Tailwind CSS, it leverages a robust Finite Element Method (FEM) engine to compute and visualize moving load effects with high precision and intuitive user experience.

## Features

- **Interactive Configuration:** Dynamically add, remove, and modify beam spans and span lengths.
- **Customizable Moving Loads:** Support for both standard truck axle configurations and lane loads with dynamic axle spacing and load values.
- **VBA-Based FEM Analysis:** Continuous Euler-Bernoulli beam analysis using the equations and load rules in [LL Analysis VBA Code.txt](LL%20Analysis%20VBA%20Code.txt), with a cached banded factorization and analysis in an inline web worker.
- **Advanced Visualizations:** Smooth, interactive SVG-based envelope charts indicating both maximum and minimum envelopes for shear, moment, and deflection.
- **Responsive Beam Schematics:** Visual representation of beam configurations and support reactions.
- **Complete Load-Case Results:** Truck, lane, and combined envelopes, upward-positive reaction summaries with continuously optimized maxima and uplift minima, and governing truck positions. Sampled reaction histories are retained in Excel, not displayed as charts.
- **Verified Lane UDL Tracer:** Automatically traces moment at the node nearest the bridge midpoint, showing favorable influence zones and partial-element intervals. Independently reassembles and solves those UDL placements to verify their demands.
- **Full-Precision Excel Export:** Input settings, span and axle configurations, all calculated cases and reaction histories, plus UDL tracer verification, intervals and influence ordinates for Lane/Envelope analyses. Deflections are not rounded to zero.

## Analysis Rules and Accuracy

The engine follows the calculation rules in the updated supplied VBA reference:

- Pinned vertical supports at every span boundary, continuous rotations, constant material properties, and cubic Hermite consistent point-load vectors.
- Element force recovery uses `k*u - f_element`; shear has two samples per element to preserve support jumps. Moments are nodal and sagging-positive; deflections are in meters; reactions are upward-positive kN.
- The truck runs in both axle orientations. For each response and placement, favorable axles are selected by the sign of their influence contribution. Automatic DLA is 40% for one selected axle, 30% for two or the original front-three group, and 25% for other groups of three or more. A multiplier `d` in [0,1] scales only DLA; a nonnegative manual override is also multiplied by `d`. Neither applies to the lane case.
- The lane case is **0.8 times the selected truck plus configurable influence-zone UDL, with no DLA** on either portion. Blank UDL defaults to 9 kN/m; zero removes only UDL, leaving the 80% lane truck. Each response has its own positive/negative loaded zones, including partial elements, and UDL may overlap the truck.
- Truck extrema are found by continuous piecewise-cubic optimization, including stationary points and one-sided element boundaries, rather than relying on the sweep grid. Both support maxima and minima/uplift use these continuous extrema; reaction-history exports remain sampled.
- Envelope mode retains separate Truck and Lane results as well as their pointwise combined envelope; use **Display load case** to inspect each.
- The automatic UDL tracer selects the moment node nearest the overall bridge midpoint (earlier node in an exact tie), **not** the governing moment location. It shows UDL-only contributions, favorable intervals and unit influence ordinates. Partial-element consistent load vectors are assembled and solved separately; a failed reconstruction raises an explicit analysis error using the VBA tolerance of `1e-9 + 1e-7 * max(abs(max), abs(min))` kNm. With zero UDL, demand and loaded intervals are zero/empty, but unit influence zones remain visible.

The stiffness matrix is factorized once with banded LDL^T (equivalent to the reference's LU). Four unit point-load solves within each element recover its exact cubic influence functions, shared across cases. Their roots partition the positive/negative UDL zones for exact integration; selected truck responses are optimized continuously in both orientations. The UI reports actual unit-load solves and element/response UDL integrations, not fictitious whole-span UDL solves. Two additional factored-system solves verify the automatic lane tracer.

### Sweep resolution and intentional safeguards

The adaptive increment considers the base step, shortest span / 40, shortest element, and smallest positive axle gap / 8. A 0.02 m minimum and a 6,000-uniform-interval cap per direction bound the solve count. The step is rounded up to 5 mm without exceeding the base step unless a safety floor or cap requires it. The effective step and its reason are displayed and exported.

Compared with the reference, the app deliberately:

- **Enforces** the step floor and count cap even when the base step is smaller (the VBA's final base-step clamp can undo these safeguards).
- Includes the end of the sweep and both signs of the reference's support-target offsets. A common exact-coordinate grid is evaluated in both axle orientations. Additional alignment samples may improve extrema slightly compared with the reference's separately sampled passes.
- Stores targeted reaction samples at their **actual lead positions**, instead of binning them into a nearby uniform position. Continuously optimized governing positions may lie between sampled history ordinates.
- Rejects invalid inputs explicitly: fractional element counts, non-finite/nonpositive stiffness and geometry, negative loads/spacings, UDL or DLA, invalid `d`, and unsupported cases. Limits are 20 axles, 12 spans for lane/envelope cases, 1,000 total elements, and 2,000,000 reaction-history ordinates per case. The web app retains support for 12 lane/envelope spans versus the VBA worksheet's 10 input rows.

Regression tests compare continuous extrema against the legacy dense LU/whole-span-pattern baseline, exact simply supported solutions and exact two-span uplift. Traced partial-element placements are also checked with an independent dense LU solver using Gauss quadrature for load-vector integration. These verify numerical implementation, not independent MIDAS Civil certification. Matching a specific MIDAS model still requires identical properties, load factors, influence-zone loading, signs, output stations and moving-load rules, followed by a project-specific benchmark. For design use, check mesh convergence and have results reviewed by a qualified engineer.

## How to Use

### Excel VBA truck selection

In [LL Analysis VBA Code.txt](LL%20Analysis%20VBA%20Code.txt), the generated **LL Input!B8** dropdown selects **CL-625** (default) or **BCL-625**. Existing input sheets gain this selector when Analyze is run. CL-625 continues to use the editable axle table unchanged; switching to BCL-625 does not overwrite that table.

BCL-625 uses axle loads **50, 140, 140, 175, 120 kN** and successive gaps **3.6, 1.2, V, 6.6 m**, as in the supplied figure. The algorithm analyzes **24 configurations**: V = **6.6, 7.1, ... , 17.6 m**, followed by **18.0 m** (a final 0.4 m step). Each configuration is optimized continuously in truck position, in both orientations. Shear, moment, deflection, optimized support maxima/minima and sampled reaction histories are enveloped over every configuration. Truck, Lane and combined Envelope cases retain their existing DLA and lane-load rules; the UDL influence zones and stiffness factorization are reused.

Results identify BCL-625, and governing truck moment/shear diagnostics report the controlling axle 3-4 spacing and reconstruct that configuration for the FEM check. **B33 remains the truck-position sampling increment**, not the variable-spacing increment. This is a discrete 0.5 m spacing search, not continuous optimization in V; check mesh/spacing convergence and obtain engineering review for design use. These options apply to the Excel VBA only; the web application is unchanged.

On Windows with Excel installed and **Trust access to the VBA project object model** enabled, run `powershell -NoProfile -File .\tests\vba-truck-spacing.test.ps1` for the Excel/VBA regression checks. The runner creates and closes its own unsaved workbook; it does not change Excel security settings.

1. **Configuration:** 
   - Set the structural material properties such as Young's Modulus ($E$) and Moment of Inertia ($I$).
   - Choose Truck, Lane, or Envelope and adjust mesh and base sweep increment.
2. **Span Setup:** 
   - Add spans and adjust their lengths in meters. The schematic will update in real time.
3. **Axle Setup:** 
   - Define custom axles and their respective loads (kN) and spacings (m).
4. **Analysis & Results:** 
   - Click **Run Analysis**. Progress is shown while the UI remains responsive. Inspect shear, moment, deflection, support summaries and the automatic lane UDL tracer, then export all calculated cases and sampled reaction histories.
   - Results and exports retain the analyzed input snapshot; changing configuration does not relabel previous results.

## Development & Running Locally

This project uses `npm` and `vite`. Use Node.js 22.18+ (tested with 24.18); the dependency-free numerical tests use Node's built-in TypeScript support. In Windows PowerShell, use `npm.cmd` if script execution policy blocks `npm`.

The app is a single npm project directly in the repository root, with one `package.json` and one `package-lock.json`. There is no nested application or npm workspace. This README is the canonical project documentation. The standalone HTML runs without installed dependencies. If dependency folders have been cleaned, restore them with `npm ci` from the root before development, linting, or rebuilding. Firebase hosting caches are generated and ignored; deployment configuration is retained.

### Project layout

- `index.html`: ready-to-open standalone app; `start.bat` opens it on Windows.
- `src/`: React UI, analysis engine and inline worker.
- `tests/`: numerical regression tests, independent VBA-equation reference, and root-layout/build-failure checks.
- `index.template.html`, `build-standalone.mjs`, `vite.config.ts`, and `tsconfig*.json`: development and build inputs.
- `dist/index.html`: generated standalone output for deployment; `firebase.json`, `.firebaserc`, and `apphosting.yaml` configure deployment from the root.

Run all npm commands from this folder, not a subfolder.

Editing `LL Analysis VBA Code.txt` does not automatically update the web solver. The web engine has been compared with the updated reference, including selected-axle DLA, continuous truck optimization, partial-element UDL zones, blank/zero UDL handling and automatic tracer verification. Excel worksheet layout/migration and VBA-only truck diagnostics are not reproduced in the web UI; further reference edits still need separate review and validation.

1. Install dependencies:
   ```bash
   npm ci
   ```
2. Start the development server:
   ```bash
   npm run dev
   ```
3. To build for production:
   ```bash
   npm run build
   ```
   This type-checks the app and rebuilds the standalone HTML at the repository root and in `dist/index.html`. Open the root `index.html` directly in a browser; the analysis worker is embedded and works without a server or internet connection. Excel export retains the existing remotely loaded SheetJS dependency and requires internet access. Vite reads `index.template.html` for development and builds without replacing the ready-to-open standalone app during development.
4. Run numerical regression tests and lint:
   ```bash
   npm test
   npm run lint
   ```
5. Serve the production build locally:
   ```bash
   npm start
   ```
   This serves `dist` on port 8080. It is separate from `start.bat`, which opens the standalone app without a server.

## Contributing

Contributions, bug reports, and feature requests are welcome! Feel free to modify, upgrade, and fork the application for educational and engineering purposes. Please refer to the [LICENSE](LICENSE) file for more information on user restrictions (e.g., selling the software itself is prohibited, but using it for paid consulting services is allowed).

## License

This project is licensed under a Custom License. See the [LICENSE](LICENSE) file for details.
