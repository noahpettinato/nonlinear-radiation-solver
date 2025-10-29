% ========================================================================
% PHASE 1 — Newton’s Method (Algorithm 1) for Pollock’s Equation (Eq. 20)
% ========================================================================
% File: newton_eq20.m
% Author: Noah Pettinato
%
% Summary:
%   Implements Newton’s method with optional damping (β ≤ 1) to solve the
%   nonlinear benchmark system f(x) = 0 defined in f_and_J.m. Serves as the
%   core iteration routine for analyzing convergence of Pollock’s Equation
%   as ε → 0⁺. Terminates on ‖f(x)‖₂ < tol or when the Jacobian becomes
%   ill-conditioned (cond(J) > 1e12).
%
% Key features:
%   • Supports damping parameter β for convergence control
%   • Uses analytic Jacobian from f_and_J.m
%   • Tracks iteration count and convergence success flag
%   • Stops if ‖f(x)‖₂ < tol or Jacobian nearly singular
%
% Inputs:
%   x0        - initial guess [x₁; x₂]
%   epsilon   - scalar parameter in f(x)
%   beta      - damping factor (0 < β ≤ 1)
%   tol       - convergence tolerance
%   max_iter  - maximum number of iterations
%
% Outputs:
%   x         - final solution approximation
%   num_iter  - number of iterations performed
%   converged - true if convergence succeeded
%
% Dependencies:
%   f_and_J.m  % returns [f, J] for given x and ε
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> [x, num_iter, converged] = newton_eq20([1; 3], 1e-6, 1, 1e-10, 100)
%
% ========================================================================

function [x, num_iter, converged] = newton_eq20(x0, epsilon, beta, tol, max_iter)
% newton_eq20   Applies Newton's method (Algorithm 1) to solve f(x) = 0.

    %% === Initialization ===
    x = x0;
    converged = false;

    %% === Newton iteration loop ===
    for k = 1:max_iter
        [f, J] = f_and_J(x, epsilon);

        if norm(f, 2) < tol
            converged = true;
            break;
        end

        if cond(J) > 1e12
            num_iter = max_iter + 1;
            return
        end

        s = J \ f;
        x = x - beta * s;
    end

    num_iter = k;
end
