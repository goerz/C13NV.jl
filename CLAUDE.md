# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

C13NV.jl is a Julia package for modeling nitrogen-vacancy (NV) centers in diamonds, focused on quantum control applications. It implements multi-level quantum systems with ground (G), excited (E), and metastable singlet (S) states, including nuclear spin coupling, magnetic field effects, and integration with the QuantumControl.jl ecosystem.

## Development Commands

### Core Development

- `make test` - Run full test suite
- `make coverage` - Run tests with coverage reporting
- `make devrepl` - Start interactive REPL with test dependencies loaded
- `make codestyle` - Format code using JuliaFormatter
- `make docs` - Build documentation

### Testing

- `julia --project=test -e 'using TestItemRunner; @run_package_tests'` - Run specific test files
- Tests are organized per-file using TestItemRunner.jl framework

### Dependencies

- Requires Julia 1.11+
- Main dependencies: QuantumControl.jl, QuantumPropagators.jl

## Architecture

### Core Modules (`src/`)

- **`models.jl`**: Contains `make_nv_system()` function that constructs quantum Hamiltonians, Liouvillians, and collapse operators for NV center systems
- **`amplitudes.jl`**: Control field parameterizations (`LinearChirp`, `ConstantDrive`) compatible with QuantumControl.jl optimization
- **`defaults.jl`**: Physical parameters for NV center systems
- **`units.jl`**: Unit system (MHz, GHz, μs, Gauss, etc.)

### Key Functions

- `make_nv_system()`: Primary function in `models.jl` for constructing quantum system representations
- State labeling uses tuples: `(N_label, ms_label, mI_label)` for readability
- System supports both Hamiltonian (unitary) and Liouvillian (open system) evolution

### Physics Conventions

- Frequencies in MHz/GHz, time in μs/ns/ms, magnetic fields in Gauss
- Complex matrices use ComplexF64 precision
- Quantum states constructed with `ket()` and `bra()` functions following standard bra-ket notation

### Integration Points

- Built for QuantumControl.jl optimization workflows
- ComponentArrays.jl for structured parameter handling in optimization problems
- QuantumPropagators.jl for time evolution computations

## Development Environment Setup

The project uses separate environments for different purposes:

- **Test environment** (`test/Project.toml`): Additional testing dependencies
- **Documentation** (`docs/Project.toml`): Documentation building tools

The `test` environment encompasses the `docs` environment. Run `make test/Manifest.toml` once to make sure the test environment is properly instantiated. Then, `julia --project=test -e …` can run Julia code. For example, to apply code formatting:


	julia --project=test -e 'using JuliaFormatter; format(["src", "docs", "test"])'

This code formatting MUST be run every time any `.jl` file is modified.

Make sure to only use explicit imports in Julia code, and that there are no imported functions or constants that are not actually used.
