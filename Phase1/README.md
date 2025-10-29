# Phase 1 — Newton’s Method Convergence for Pollock’s Equation

## Objective
Investigate the convergence behavior of Newton’s method when applied to Pollock’s nonlinear benchmark system (Eq. 20).  
This phase explores how the convergence domains (basins of attraction) evolve as the small parameter ε → 0⁺ and validates the effect of damping and Jacobian conditioning on convergence.

---

## Governing System (Pollock’s Equation, Eq. 20)
Let \( x = [x_1; x_2] \) and \( y = x_2 - 3 \):
\begin{cases}
f_1(x) = (x_1 - 1) + y^2, \\
f_2(x) = \varepsilon y + \tfrac{3}{2}(x_1 - 1)y + y^2 + y^3.
\end{cases}
Newton’s method iterates
x_{k+1} = x_k - \beta J_f(x_k)^{-1} f(x_k),
where \( \beta \le 1 \) is a damping parameter.

---

## Files Included
| File | Description |
|------|--------------|
| **`f_and_J.m`** | Defines nonlinear system *f(x)* and analytic Jacobian *J_f(x)* for Pollock’s Equation. |
| **`newton_eq20.m`** | Core Newton solver with optional damping (β ≤ 1) and Jacobian conditioning check. |
| **`convergence_plot.m`** | Visualizes Newton convergence domain for a fixed ε = 10⁻⁶. |
| **`compare_convergence_epsilons.m`** | Compares convergence domains across several fixed ε values (10⁻³–10⁻⁶). |
| **`animate_convergence_epsilons.m`** | Animates how the convergence basin evolves as ε decreases logarithmically (10⁻¹ → 10⁻⁶). |

---

## Numerical Methods
- **Root-Finding Algorithm:** Newton’s method with optional damping.  
- **Convergence Test:** \( \|F(x)\|_\infty < 10^{-10} \), max 100 iterations.  
- **Visualization:**
  - Iteration count mapped over a 2D grid of initial guesses (x₁,x₂).
  - Color-coded heatmaps and animated evolution as ε → 0⁺.
- **Equilibria:** Two steady states (x₊, x₋) appear as ε → 0⁺.

---

## Expected Outputs
- Static convergence domain plots for several ε values.  
- Animated sequence showing deformation of the convergence basins as ε → 0⁺.  
- Iteration count color map highlighting divergence regions.

---

## How to Run
In MATLAB:
```matlab
>> convergence_plot
>> compare_convergence_epsilons
>> animate_convergence_epsilons
