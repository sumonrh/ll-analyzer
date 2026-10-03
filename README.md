# LL Analyzer

LL Analyzer is a modern, professional web application designed for comprehensive Live Load Analysis of continuous beams and out-of-the-box NU Girder designs. Built with React, Vite, and Tailwind CSS, it leverages a robust Finite Element Method (FEM) engine to compute and visualize moving load effects with high precision and intuitive user experience.

## Features

- **Interactive Configuration:** Dynamically add, remove, and modify beam spans and span lengths.
- **Customizable Moving Loads:** Support for both standard truck axle configurations and lane loads with dynamic axle spacing and load values.
- **VBA-Based FEM Analysis:** Continuous Euler-Bernoulli beam analysis using the equations and load rules in [LL Analysis VBA Code.txt](LL%20Analysis%20VBA%20Code.txt), with a cached banded factorization and analysis in an inline web worker.
- **Advanced Visualizations:** Smooth, interactive SVG-based envelope charts indicating both maximum and minimum envelopes for shear, moment, and deflection.
- **Responsive Beam Schematics:** Visual representation of beam configurations and support reactions.
- **Complete Load-Case Results:** Truck, lane, and combined envelopes, upward-positive reaction summaries, exact governing truck positions, individual support reaction diagrams, and an all-supports diagram.
- **Full-Precision Excel Export:** Input settings, span and axle configurations, and all calculated cases, including reaction histories. Deflections are not rounded to zero.

## Analysis Rules and Accuracy

The engine follows the supplied VBA reference:

- Pinned vertical supports at every span boundary, continuous rotations, constant material properties, and cubic Hermite consistent point-load vectors.
- Element force recovery uses `k*u - f_element`; shear has two samples per element to preserve support jumps. Moments are nodal and sagging-positive; deflections are in meters; reactions are upward-positive kN.
- The truck runs in both axle orientations. Automatic DLA is evaluated on **every span**, governing with 40% for one axle, 30% for two, and 25% for three or more. A nonnegative user override applies only to the truck-only case.
- The lane case is **0.8 times the truck plus 9 kN/m patterned span UDL, with no DLA** on either portion. Loaded/unloaded patterns include the empty pattern. Each response ordinate has its own governing pattern.
- Envelope mode retains separate Truck and Lane results as well as their pointwise combined envelope; use **Display load case** to inspect each.

The app uses mathematically equivalent optimizations rather than repeating unnecessary solves: the stiffness matrix is factorized once with banded LDL^T (equivalent to the reference's LU); the unscaled truck sweep is shared between cases; and positive/negative contributions from individual span UDL responses give the exact extrema of all `2^n` patterns in just `n` solves.

### Sweep resolution and intentional safeguards

The adaptive increment considers the base step, shortest span / 40, shortest element, and smallest positive axle gap / 8. A 0.02 m minimum and a 6,000-uniform-interval cap per direction bound the solve count. The step is rounded up to 5 mm without exceeding the base step unless a safety floor or cap requires it. The effective step and its reason are displayed and exported.

Compared with the reference, the app deliberately:

- **Enforces** the step floor and count cap even when the base step is smaller (the VBA's final base-step clamp can undo these safeguards).
- Includes the end of the sweep and both signs of the reference's support-target offsets. A common exact-coordinate grid is evaluated in both axle orientations. Additional alignment samples may improve extrema slightly compared with the reference's separately sampled passes.
- Stores targeted reaction samples at their **actual lead positions**, instead of binning them into a nearby uniform position; the governing position and diagram ordinate remain consistent.
- Rejects invalid inputs explicitly: fractional element counts, non-finite/nonpositive stiffness and geometry, negative loads/spacings or DLA, and unsupported cases. Limits are 20 axles, 12 spans for lane/envelope cases, 10,000 total elements, and 2,000,000 reaction-history ordinates per case.

Regression tests compare the optimized engine against an independent dense LU implementation of the VBA equations, exhaustive lane patterns, and exact simply supported point-load/UDL solutions. These verify numerical implementation, not independent MIDAS Civil certification. Matching a specific MIDAS model still requires identical properties, load factors, lane patterning, signs, output stations and moving-load resolution, followed by a project-specific benchmark. For design use, check mesh/sweep convergence and have results reviewed by a qualified engineer.

## How to Use

1. **Configuration:** 
   - Set the structural material properties such as Young's Modulus ($E$) and Moment of Inertia ($I$).
   - Choose Truck, Lane, or Envelope and adjust mesh and base sweep increment.
2. **Span Setup:** 
   - Add spans and adjust their lengths in meters. The schematic will update in real time.
3. **Axle Setup:** 
   - Define custom axles and their respective loads (kN) and spacings (m).
4. **Analysis & Results:** 
   - Click **Run Analysis**. Progress is shown while the UI remains responsive. Inspect shear, moment, deflection, support summaries and reaction diagrams, then export all calculated cases.
   - Results and exports retain the analyzed input snapshot; changing configuration does not relabel previous results.

## Development & Running Locally

This project uses `npm` and `vite`. Use Node.js 22.18+ (tested with 24.18); the dependency-free numerical tests use Node's built-in TypeScript support. In Windows PowerShell, use `npm.cmd` if script execution policy blocks `npm`.

1. Install dependencies:
   ```bash
   npm install
   ```
2. Start the development server:
   ```bash
   npm run dev
   ```
3. To build for production:
   ```bash
   npm run build
   ```
   This type-checks the app and rebuilds the standalone HTML at the repository root and in `ll-analyzer`, plus the deployable `ll-analyzer/dist` output. Open either standalone `index.html` directly in a browser; the analysis worker is embedded and works without a server or internet connection. Excel export retains the existing remotely loaded SheetJS dependency and requires internet access.
4. Run numerical regression tests and lint:
   ```bash
   npm test
   npm run lint
   ```

## Contributing

Contributions, bug reports, and feature requests are welcome! Feel free to modify, upgrade, and fork the application for educational and engineering purposes. Please refer to the [LICENSE](LICENSE) file for more information on user restrictions (e.g., selling the software itself is prohibited, but using it for paid consulting services is allowed).

## License

This project is licensed under a Custom License. See the [LICENSE](LICENSE) file for details.
