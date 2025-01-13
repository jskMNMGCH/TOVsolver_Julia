# TOVsolver_Julia
This project contains Julia scripts to solve the Tolman–Oppenheimer–Volkoff (TOV) equations, which are fundamental in modeling the structure of compact stars such as neutron stars. The solver calculates key parameters like mass, radius, and tidal deformability based on input energy density and pressure profiles.

## Project Structure

### Files

- **`main.jl`**: The main module that provides high-level functions for processing input data and running the solver.
- **`solver_code.jl`**: A module containing core functions for solving the TOV equations, debugging, and calculating tidal deformabilities.
- **`constants.jl`**: (Not provided here) Assumed to define physical and unit conversion constants used in the calculations.

### Modules

- **`MainModule`**: Provides functions for preprocessing data (e.g., ensuring monotonicity) and orchestrates the execution of the TOV solver.
- **`SolverCode`**: Implements the TOV equations, numerical solvers, and utility functions for debugging and tidal deformability calculations.

## Usage

### Setup
1. Install Julia from [the official website](https://julialang.org/).
2. Ensure the following Julia packages are installed:
   - `DifferentialEquations`
   - `Dates`

   Install packages using:
   ```julia
   using Pkg
   Pkg.add("DifferentialEquations")
   Pkg.add("Dates")
   ```

3. Place all project files in the same directory to ensure correct module inclusion.

### Running the Solver
1. Open the Julia REPL or a Julia-supported IDE (e.g., Juno or VS Code).
2. Load the main module:
   ```julia
   include("main.jl")
   using .MainModule
   ```

3. Prepare input data:
   - Energy densities (`ε`) and pressures (`P`) as arrays.

4. Call the main functions:
   - **`make_eos_monotonic(e, P)`**: Ensures monotonicity of the equation of state.
   - **`out_RMT(ε, pres; debug=false)`**: Runs the TOV solver and returns radius, mass, and tidal deformability.

   Example:
   ```julia
   ε = [1e-5, 1e-4, 1e-3, ...]  # Example energy densities
   P = [1e-6, 1e-5, 1e-4, ...]  # Example pressures

   mono_e, mono_P = MainModule.make_eos_monotonic(ε, P)
   results, solutions = MainModule.out_RMT(mono_e, mono_P, debug=true)
   ```

### Debugging
- **Debugging Log**: Logs detailed calculations to a timestamped file in the `debug/` directory. The debug data includes interpolated energy density, pressure, and other parameters for troubleshooting.
- Enable debugging by setting `debug=true` in `out_RMT`.

## Key Functions

### In `main.jl`

#### `make_eos_monotonic(e, P)`
- **Description**: Removes non-monotonic points from energy density (`e`) and pressure (`P`) arrays.
- **Arguments**:
  - `e`: Array of energy densities.
  - `P`: Array of pressures.
- **Returns**: Monotonic arrays of energy densities and pressures.

#### `out_RMT(ε, pres; debug=false)`
- **Description**: Solves the TOV equations for given energy densities (`ε`) and pressures (`pres`).
- **Arguments**:
  - `ε`: Array of energy densities.
  - `pres`: Array of pressures.
  - `debug`: Enables detailed logging if set to `true`.
- **Returns**: A tuple containing:
  - Radii (`R`), masses (`M`), and tidal deformabilities (`Λ`).
  - Solution objects for detailed inspection.

### In `solver_code.jl`

#### `solveTOV_RMT(center_idx, ε, pres, debug_flag)`
- **Description**: Numerically solves the TOV equations for a given central density and pressure.
- **Arguments**:
  - `center_idx`: Index of the central density in input arrays.
  - `ε`: Array of energy densities.
  - `pres`: Array of pressures.
  - `debug_flag`: Enables logging if set to `true`.
- **Returns**: A tuple containing:
  - Radius (`R`), mass (`M`), and tidal deformability (`Λ`).
  - Detailed solution object.

#### `tidal_deformability(y, M, R)`
- **Description**: Calculates tidal deformability for a given compactness and other parameters.

#### `TOV_def!(du, u, p, t)`
- **Description**: Defines the system of differential equations for the TOV equations.

## Outputs
- **Radii (`R`)**: In kilometers.
- **Masses (`M`)**: In solar mass units.
- **Tidal Deformabilities (`Λ`)**: Dimensionless.
- **Debug Logs**: Saved in the `debug/` directory for troubleshooting.

## Example Workflow
1. Define the energy density (`ε`) and pressure (`P`) profiles.
2. Ensure monotonicity using `make_eos_monotonic`.
3. Call `out_RMT` to solve the TOV equations and retrieve the results.
4. Analyze the output radii, masses, and tidal deformabilities for your specific application.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments
This solver leverages the Julia `DifferentialEquations` package for numerical integration and the `Dates` package for logging. It is designed for research applications in astrophysics and compact star modeling.

