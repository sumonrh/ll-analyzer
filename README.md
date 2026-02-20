# LL Analyzer

LL Analyzer is a modern, professional web application designed for comprehensive Live Load Analysis of continuous beams and out-of-the-box NU Girder designs. Built with React, Vite, and Tailwind CSS, it leverages a robust Finite Element Method (FEM) engine to compute and visualize moving load effects with high precision and intuitive user experience.

## Features

- **Interactive Configuration:** Dynamically add, remove, and modify beam spans and span lengths.
- **Customizable Moving Loads:** Support for both standard truck axle configurations and lane loads with dynamic axle spacing and load values.
- **Real-time FEM Analysis:** Instant calculation of bending moments, shear forces, deflections, and support reactions as load conditions change.
- **Advanced Visualizations:** Smooth, interactive SVG-based envelope charts indicating both maximum and minimum envelopes for shear, moment, and deflection.
- **Responsive Beam Schematics:** Visual representation of beam configurations and support reactions.

## How to Use

1. **Configuration:** 
   - Set the structural material properties such as Young's Modulus ($E$) and Moment of Inertia ($I$).
   - Choose the load case (Truck or Lane load) and adjust analysis increments.
2. **Span Setup:** 
   - Add spans and adjust their lengths in meters. The schematic will update in real time.
3. **Axle Setup:** 
   - Define custom axles and their respective loads (kN) and spacings (m).
4. **Analysis & Results:** 
   - Click the **Run Analysis** button. The application will solve the continuous beam system and present comprehensive envelope charts for shear, moment, and deflection, along with maximum and minimum reaction diagrams.

## Development & Running Locally

This project uses `npm` as the package manager and `vite` for fast local development.

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

## Contributing

Contributions, bug reports, and feature requests are welcome! Feel free to modify, upgrade, and fork the application for educational and non-commercial purposes. Please refer to the [LICENSE](LICENSE) file for more information on user restrictions.

## License

This project is licensed under a Custom Non-Commercial License. See the [LICENSE](LICENSE) file for details.
