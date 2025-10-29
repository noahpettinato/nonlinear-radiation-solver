# Phase 3 — Time-Dependent Diffusion with Time-Lagged Nonlinear Robin Boundary

## Objective
Simulate the transient heat equation with a **nonlinear radiation-type Robin boundary condition**, using fully implicit (Backward Euler) time-stepping.  
This phase introduces temporal evolution and compares two spatial boundary treatments:
1. **FD (Backward Difference)** — first-order accurate,
2. **CCFD (Ghost Point)** — second-order accurate.

The radiation term is evaluated using a **time-lagged** boundary value from the previous timestep, simplifying the nonlinear flux evaluation.

---

## Governing PDE

$$
c\,u_t = k\,u_{xx}, \quad x \in (0,1),
$$

with boundary conditions

$$
u(0,t) = 0, \qquad
-\,k\,u_x(1,t) =  lpha\,[u(1,t-\Delta t)^4 - u_*^4] + g.
$$

Here \(c\) is the heat capacity, \(k\) the conductivity, \( lpha\) the radiation coefficient, \(u_*\) the reference equilibrium temperature, and \(g\) the source term.

---

## Files Included
| File | Description |
|------|--------------|
| **`phase3_fd_step.m`** | Advances one time step using **Backward Euler** in time and **backward difference** at \(x=1\) (first-order BC). |
| **`phase3_ccfd_step.m`** | Advances one time step using **Backward Euler** in time and **ghost-point centered scheme** at \(x=1\) (second-order BC). |
| **`robin_phase3.m`** | Main driver script. Evolves both FD and CCFD solutions in time and compares results at \(T_{	ext{final}}\). |

---

## Numerical Methods
- **Temporal Discretization:** Backward Euler (first-order implicit).  
- **Spatial Discretizations:**
  - FD (one-sided) — first-order enforcement of Robin flux.  
  - CCFD (ghost point) — second-order enforcement of Robin flux.  
- **Nonlinearity Handling:** Radiation term evaluated explicitly using \(u(1,t-\Delta t)\) (time-lagged).  
- **Initial Condition:** \( u(x,0) = x \).  
- **Boundary Conditions:**  
  \(u(0,t) = 0,\)  \( -k u_x(1,t) =  lpha(u(1,t-\Delta t)^4 - u_*^4) + g.\)

---

## Expected Outputs
- Printed summary:
  ```
  alpha = 1.0, ||FD - CCFD||_inf = 2.13e-04
  ```
- Comparison plots:
  - Final-time profiles of \(u(x,T_{	ext{final}})\) for FD vs CCFD.
  - Clear visual demonstration that the CCFD scheme maintains higher accuracy.

---

## How to Run
In MATLAB:
```matlab
>> robin_phase3
```

The script automatically runs for multiple α-values (e.g., 0.1, 1.0, 10.0) and generates plots comparing FD and CCFD methods.

---

## Dependencies
| Function | Used By |
|-----------|----------|
| `phase3_fd_step.m` | `robin_phase3.m` |
| `phase3_ccfd_step.m` | `robin_phase3.m` |
| (Shared functions from Phase 2 may be reused for consistency.) | — |

No external toolboxes are required.

---

## Notes
Phase 3 extends the steady-state framework (Phase 2) to the **time-dependent** case.  
By comparing FD and CCFD formulations under implicit time integration, this phase confirms that the ghost-point method provides improved spatial accuracy while maintaining numerical stability in transient simulations.
