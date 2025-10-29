# Phase 2 — Nonlinear Steady-State Robin Boundary Problem

## Objective
Develop and compare finite difference discretizations for the **nonlinear radiation-type Robin boundary condition** that arises in steady-state heat transfer.  
This phase extends the linear benchmark (Phase 0) by introducing a nonlinear flux term at the boundary and assembling analytic Jacobians for Newton’s method.

---

## Governing Equation

$$
-\,k\,u_{xx} = f(x), \quad x \in (0,1)
$$

with boundary conditions

$$
u(0) = 0, \qquad
-\,k\,u_x(1) = \alpha\,(u(1)^4 - u_*^4) + g.
$$

Here \( \alpha \) is the radiation coefficient, \( g \) a boundary source term, and \( u_* \) a prescribed equilibrium temperature.

---

## Files Included
| File | Description |
|------|--------------|
| **`f_and_J_robin_fd.m`** | Constructs residual **F(U)** and analytic Jacobian **JF(U)** using a first-order backward-difference (FD) discretization at \(x=1\). |
| **`f_and_J_robin_ccfd.m`** | Constructs **F(U)** and **JF(U)** using the second-order **centered control finite difference (CCFD)** scheme with a ghost point. |
| **`newton_eq20.m`** | General Newton solver with adaptive damping and analytic Jacobian input. |
| **`robin_phase2.m`** | Runs steady-state Newton solves for both FD and CCFD formulations at several values of \( \alpha \). |
| **`robin_phase2_convergence.m`** | Performs grid-refinement studies for each \( \alpha \), verifying first- and second-order accuracy in \( \|U - u_{\text{ref}}\|_\infty \). |

---

## Numerical Methods
- **Spatial Discretizations**
  - **FD (One-Sided):** First-order backward difference for \(u_x(1)\).
  - **CCFD (Ghost Point):** Second-order centered difference enforcing the nonlinear Robin flux.
- **Nonlinearity Handling:** Solved using **Newton’s method** with analytic Jacobians.
- **Reference Solution:** Linear profile \( u_{\text{ref}}(x) = A_{\text{ref}}x \), where \(A_{\text{ref}}\) satisfies  
  \(-kA_{\text{ref}} = \alpha(A_{\text{ref}}^4 - u_*^4)\).
- **Error Metric:**  
  \[
  E_\infty = \|U - u_{\text{ref}}\|_\infty.
  \]

---

## Expected Outputs
- Printed Newton convergence history (iteration count, residual norms).  
- Convergence plots \(E_\infty\) vs \(h\) showing:
  - FD → O(h) behavior  
  - CCFD → O(h²) behavior  
- Comparison figures of FD and CCFD steady-state profiles for various \( \alpha \).

---

## How to Run
In MATLAB:
```matlab
>> robin_phase2
>> robin_phase2_convergence
