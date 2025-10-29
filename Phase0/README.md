# Phase 0 — Linear Steady-State Problem with Robin Boundary Conditions

## Objective
Establish and verify the numerical accuracy of first-order (D₋) and second-order (D₀) finite difference schemes for the 1D steady-state heat equation with linear Robin boundary conditions.

The results provide a controlled, linear benchmark before introducing nonlinear radiation effects in later phases.

---

## Governing Equation
- u''(x) = f(x), \quad x \in (0,1)
with Robin boundary conditions
u'(0) + \alpha_0 u(0) = g_0, \quad -u'(1) = \alpha_1 u(1) + g_1.

Exact solution: \( u(x) = e^{2x} \)

---

## Files Included
| File | Description |
|------|--------------|
| **`robin.m`** | Main script. Implements both D₋ (one-sided) and D₀ (ghost-point) finite difference schemes, computes errors, and plots convergence rates. |

---

## Numerical Methods
- **D₋ (One-sided):** First-order scheme using forward/backward differences at the boundaries.  
- **D₀ (Centered + Ghost Point):** Second-order scheme enforcing Robin BCs via ghost-point elimination.  
- **Error Metrics:** Infinity norm and discrete L² norm.  
- **Convergence Verification:** Expected rates of O(h) for D₋ and O(h²) for D₀.

---

## Expected Outputs
- Console table of error norms and convergence rates.  
- Two log–log plots:
  - ∞-norm error vs. grid size \( h \)
  - L²-norm error vs. grid size \( h \)

---

## How to Run
In MATLAB:
```matlab
>> robin
