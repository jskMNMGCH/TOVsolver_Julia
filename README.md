# README

## Overview

This project provides a Julia implementation for solving the Tolman–Oppenheimer–Volkoff (TOV) equations and related astrophysical calculations. It uses a piecewise polytropic equation of state (EoS) to model neutron stars, calculates their mass, radius, and tidal deformability, and includes utilities for managing the EoS data.

## Files and Modules

### 1. `main.jl`
- **Purpose**: Acts as the main entry point for calculations.
- **Key Functions**:
  - `make_eos_monotonic(e, P)`: Processes energy density (`e`) and pressure (`P`) arrays to ensure monotonicity by removing non-physical regions.
  - `out_RMT(ε, pres; ...)`: Solves the TOV equations to compute radius, mass, and tidal deformability for given energy densities (`ε`) and pressures (`pres`).
- **Dependencies**: 
  - `solver_code.jl` for solving the TOV equations.

### 2. `solver_code.jl`
- **Purpose**: Contains the numerical solvers and tools to compute stellar properties.
- **Key Functions**:
  - `Debug(input_list; base_filename)`: Logs data for debugging.
  - `tidal_deformability(y, M, R)`: Computes tidal deformability given dimensionless compactness and stellar properties.
  - `solveTOV_RMT(center_idx, ε, pres, debug_flag)`: Solves the TOV equations for a given central density and pressure using an ODE solver.
  - `TOV_def!(...)`: Defines the TOV differential equations.
- **Dependencies**: Relies on `DifferentialEquations.jl` for ODE solving.

### 3. `piecewise_polytrope_eos_cgs.jl`
- **Purpose**: Implements the piecewise polytropic EoS framework in CGS units.
- **Key Features**:
  - Defines fixed crust parameters for the EoS.
  - Provides functions to calculate pressure, energy density, and polytropic constants for different density segments.
- **Key Functions**:
  - `p_rho_def(rho, K, Gamma)`: Calculates pressure from density using the polytropic relation.
  - `eps_rho_def(rho, K, Gamma, a)`: Computes energy density for given density and model parameters.
  - `calc_outerEoS()`: Computes parameters for the outer crust.
  - `make_polyEos(rho_lim_arr, a_arr, ...)`: Constructs energy density and pressure arrays for the entire density range.

### 4. `constants.jl`
- **Purpose**: Contains physical and conversion constants used throughout the code.

## Usage

### Prerequisites
- Julia version 1.6 or higher.
- Install required Julia packages:
  ```julia
  using Pkg
  Pkg.add(["DifferentialEquations", "Dates"])
  ```

### Running the Code
1. Include the necessary modules in your Julia environment:
   ```julia
   include("main.jl")
   using .MainModule
   ```
2. Prepare energy density and pressure data.
3. Call `out_RMT` to calculate stellar properties:
   ```julia
   results, solutions = MainModule.out_RMT(ε, pres)
   ```

## Key Concepts

### 1. Tolman–Oppenheimer–Volkoff (TOV) Equations
The TOV equations describe the structure of a spherically symmetric star in hydrostatic equilibrium, accounting for general relativity.

### 2. Piecewise Polytropic EoS
The piecewise polytropic model approximates the neutron star EoS by dividing it into segments, each described by a polytropic relation:
\[ P(\rho) = K \rho^\Gamma \]

### 3. Tidal Deformability
Tidal deformability quantifies how easily a neutron star deforms under an external gravitational field, crucial for gravitational wave studies.

## References
- Read et al., "Constraints on a phenomenologically parameterized neutron-star equation of state," Phys. Rev. D, 79, 124032 (2009). DOI: [10.1103/PhysRevD.79.124032](https://doi.org/10.1103/PhysRevD.79.124032)

## Contributing
Feel free to submit issues or pull requests to improve this project.

## License
This project is licensed under the MIT License. See the LICENSE file for details.