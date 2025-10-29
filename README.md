# Nonlinear Radiation Capsule (Version 1.0)

### Author
Noah Pettinato  
Code & Documentation Editor — in collaboration with Professor Malgorzata Peszynska  

## Overview
This MATLAB capsule investigates one-dimensional heat conduction with nonlinear radiation-type Robin boundary conditions.  
The project evolves across five phases of increasing complexity—from a linear steady-state problem to a fully implicit nonlinear solver.

Spatial discretization is performed with:
- Finite Difference (FD): first-order boundary enforcement  
- Centered Control Finite Difference (CCFD): second-order ghost-point formulation  

Time integration uses the Backward Euler method.  
Nonlinear terms are handled with Newton’s method using analytic Jacobians and optional damping.

## Project Structure
```
Project/
├── README.md
├── Phase0/   # Linear steady-state Robin problem
├── Phase1/   # Pollock benchmark for Newton’s method
├── Phase2/   # Nonlinear steady-state Robin problem
├── Phase3/   # Time-lagged nonlinear Robin BC
└── Phase4/   # Fully implicit nonlinear radiation problem
```

Each phase folder contains its own MATLAB scripts, analytic Jacobian builders, and driver routines.  
See individual README_Phase*.md files for phase-specific details.

## Governing Model

### PDE
$$
c\,u_t = k\,u_{xx} + f(x), \qquad x \in (0,1),\ t>0
$$

### Boundary Conditions
$$
u(0,t) = 0, \qquad
-\,k\,u_x(1,t) =  lpha\,(u(1,t)^4 - u_*^4) + g.
$$

### Initial Condition (for transient cases)
$$
u(x,0) = u_{init}(x).
$$

## Phase Summaries

| Phase | Description | Solver | Boundary Scheme | Time Stepping |
|:------|:-------------|:--------|:----------------|:---------------|
| 0 | Linear steady-state with Robin BC | Direct | FD, CCFD | — |
| 1 | Pollock benchmark system (Newton validation) | Newton | — | — |
| 2 | Steady nonlinear radiation BC | Newton | FD, CCFD | — |
| 3 | Time-lagged nonlinear Robin BC | Linear per step | FD, CCFD | Backward Euler |
| 4 | Fully implicit nonlinear radiation PDE | Newton per step | FD, CCFD | Backward Euler |

## Numerical Methods

- Spatial Discretization
  - Interior: second-order centered finite differences.
  - Right boundary:  
    - FD (backward difference, O(h))  
    - CCFD (ghost-point centered, O(h²))
- Time Discretization: fully implicit Backward Euler.  
- Nonlinearity Handling: Newton’s method with analytic Jacobians.  
- Convergence Criterion:  
  $$ \|F(U)\|_2 < 10^{-8}. $$
- Damping: Adaptive halving strategy when residual increases.

## Implementation Highlights

- Modular design separates residual/Jacobian assembly from solver drivers.  
- newton_eq20.m reused across all nonlinear phases.  
- Boundary flux terms encapsulated in reusable radiation modules.  
- Each driver script prints iteration logs and generates comparison plots automatically.  

## Key Results

- Phase 0: verified O(h) and O(h²) convergence for FD and CCFD.  
- Phase 1: validated Newton’s stability and convergence on Pollock’s system.  
- Phase 2: accurate enforcement of nonlinear radiation BC via Newton.  
- Phase 3: stable transient evolution using lagged flux; FD and CCFD converge to same steady state.  
- Phase 4: fully implicit nonlinear solver stable even for large α and Δt; confirms asymptotic agreement with Phase 2 steady state.

## How to Run

In MATLAB:
```matlab
% Phase 0 – Linear Robin test
robin

% Phase 1 – Pollock benchmark
convergence_plot
compare_convergence_epsilons
animate_convergence_epsilons

% Phase 2 – Nonlinear steady-state
robin_phase2
robin_phase2_convergence

% Phase 3 – Time-lagged radiation BC
robin_phase3

% Phase 4 – Fully implicit nonlinear PDE
phase4_fd
phase4_ccfd
phase4_compare
```

Each script prints diagnostics and generates figures in the MATLAB workspace.

## References
1. Sara Pollock and Hunter Schwartz. *Benchmarking Results for the Newton–Anderson Method.* Results in Applied Mathematics, 8 (2020): 100095.  
2. Malgorzata Peszynska. *Modeling Nonlinear Radiation Boundary Conditions in Heat Conduction Problems.* Lecture Notes and Project Instructions, Oregon State University, 2025.  
3. Randall J. LeVeque. *Finite Difference Methods for Ordinary and Partial Differential Equations: Steady-State and Time-Dependent Problems.* SIAM, 2007.  
4. C. T. Kelley. *Solving Nonlinear Equations with Newton’s Method.* SIAM, 2003.  
5. J. David Logan. *Applied Mathematics.* 3rd Edition, Wiley, 2006.  
6. MATLAB Documentation. MathWorks, 2025. https://www.mathworks.com/help/  
7. Notes from MTH 452 and MTH 453 (Numerical ODEs and PDEs). Oregon State University, 2024–2025.

## Conclusion
This capsule provides a complete solver pipeline for 1D heat equations with nonlinear radiation effects.  
Each phase builds conceptually and computationally upon the previous one, culminating in a validated, modular MATLAB framework for nonlinear parabolic PDEs with complex boundary physics.

Version: 1.0    Date: August 2025  
Institution: Oregon State University
