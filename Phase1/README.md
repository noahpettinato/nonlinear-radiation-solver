# Phase 1 – Pollock Benchmark (Newton's Method Validation)

This phase validates the nonlinear solver implementation using Pollock’s benchmark equation (Equation 20 in the capsule), ensuring Newton’s method with analytic Jacobians converges as expected.

### Problem Overview
The system models a simplified nonlinear reaction term with known analytic steady-state. This phase is designed to verify the robustness of Newton’s method before applying it to the radiation boundary problem.

### Key Features
- Implements Newton’s method with analytic Jacobian construction.
- Demonstrates quadratic convergence near the true solution.
- Serves as a validation test for later nonlinear PDE solvers.

### Key Files
- `pollock_fd.m` – Newton solver applied to Pollock benchmark

### Outputs
- Iteration convergence plot
- Console report of Newton iteration count and residual norms

### Usage
```matlab
>> pollock_fd
