# Phase 4 — Fully Implicit Time-Dependent Nonlinear Robin Problem

## Objective
Implement a **fully implicit time-stepping scheme** for the nonlinear heat equation with radiation-type Robin boundary conditions.  
This phase advances the Phase 3 formulation by evaluating the nonlinear boundary term implicitly at the current time step and solving the coupled system using **Newton’s method** with analytic Jacobians.

---

## Governing PDE

$$
c\,u_t = k\,u_{xx}, \quad x \in (0,1),
$$

with boundary conditions

$$
u(0,t) = 0, \qquad
-\,k\,u_x(1,t) =  lpha\,(u(1,t)^4 - u_*^4) + g.
$$

Here \(c\) is the heat capacity, \(k\) the conductivity, \( lpha\) the radiation coefficient, \(u_*\) the reference equilibrium temperature, and \(g\) the boundary source term.

---

## Files Included
| File | Description |
|------|--------------|
| **`f_and_J_phase4_fd.m`** | Residual and Jacobian for fully implicit scheme using **backward difference (FD)** at \(x=1\). |
| **`f_and_J_phase4_ccfd.m`** | Residual and Jacobian for fully implicit scheme using **centered (CCFD)** method with a ghost point at \(x=1\). |
| **`newton_eq20.m`** | Damped Newton solver with adaptive step size and analytic Jacobian input. |
| **`phase4_fd.m`** | Time-dependent solver using **FD** boundary scheme. |
| **`phase4_ccfd.m`** | Time-dependent solver using **CCFD** boundary scheme. |
| **`phase4_compare.m`** | Runs both FD and CCFD methods, compares final-time solutions, and prints \( \|U_{FD} - U_{CCFD}\|_\infty \). |

---

## Numerical Methods

- **Time Discretization:** Fully implicit Backward Euler.
- **Spatial Discretizations:**
  - FD (one-sided) — first-order boundary enforcement.
  - CCFD (ghost point) — second-order boundary enforcement.
- **Nonlinearity Handling:** Newton’s method with analytic Jacobian from `f_and_J_phase4_fd.m` or `f_and_J_phase4_ccfd.m`.
- **Convergence Tolerance:**  
  $$ \|F(U)\|_2 < 10^{-8}. $$
- **Initial Condition:**  
  $$ u(x,0) = \sin(\pi x). $$

---

## Expected Outputs
- Printed Newton iteration logs (residual norms per timestep).  
- Warnings for non-convergent steps (Jacobian conditioning checks).  
- Final-time comparison plot between FD and CCFD schemes.  
- Reported infinity-norm error:  
  ```
  || FD - CCFD ||_inf = 3.42e-05
  ```

---

## How to Run
In MATLAB:
```matlab
>> phase4_fd
>> phase4_ccfd
>> phase4_compare
```
Each script runs a fully implicit time integration with nonlinear boundary conditions.

---

## Dependencies
| Function | Used By |
|-----------|----------|
| `f_and_J_phase4_fd.m` | `phase4_fd.m`, `phase4_compare.m` |
| `f_and_J_phase4_ccfd.m` | `phase4_ccfd.m`, `phase4_compare.m` |
| `newton_eq20.m` | All solver scripts |

No external toolboxes are required.

---

## Notes
Phase 4 represents the **final and most complete implementation** of the nonlinear radiation problem.  
The FD and CCFD schemes are now compared under a **fully coupled implicit formulation**, demonstrating both stability and second-order spatial accuracy for the CCFD approach.
