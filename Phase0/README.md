# Phase 0 – Linear Steady-State Problem with Robin Boundary Condition

This phase establishes a baseline accuracy comparison for the steady-state linear problem

- k u_{xx} = 0, \quad u(0) = 0, \quad -k u_x(1) = \alpha (u(1) - u_*).

Two spatial discretizations are implemented:
1. **FD (Finite Difference)** – backward one-sided derivative at \(x=1\)
2. **CCFD (Centered Control Finite Difference)** – ghost-point formulation for second-order accuracy

### Key Features
- Linear test problem used to verify implementation of boundary treatment.
- Confirms second-order accuracy of the CCFD method.
- Serves as the foundation for later nonlinear extensions.

### Key Files
- `robin_fd.m` – one-sided (FD) method
- `robin_ccfd.m` – ghost-point (CCFD) method

### Outputs
- Printed infinity-norm error between FD and CCFD
- Plot of \(u(x)\) for both methods

### Usage
```matlab
>> robin_fd
>> robin_ccfd
