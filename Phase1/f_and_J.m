% ========================================================================
% PHASE 1 — Nonlinear System and Jacobian for Pollock’s Equation (Eq. 20)
% ========================================================================
% File: f_and_J.m
% Author: Noah Pettinato
%
% Summary:
%   Defines the nonlinear benchmark system f(x) = 0 and its Jacobian J_f(x)
%   used to test Newton’s method convergence for Pollock’s Equation. The
%   parameter ε controls nonlinearity strength, with ε → 0⁺ producing a
%   near-singular Jacobian at the equilibrium roots (x₊, x₋).
%
% Mathematical Form:
%   Let x = [x₁; x₂], and y = x₂ − 3:
%       f₁(x) = (x₁ − 1) + y²
%       f₂(x) = εy + (3/2)(x₁ − 1)y + y² + y³
%
% Key outputs:
%   • f — vector-valued residual function f(x)
%   • J — 2×2 Jacobian matrix of f(x)
%
% Inputs:
%   x        - 2×1 vector [x₁; x₂]
%   epsilon  - scalar parameter controlling nonlinearity
%
% Outputs:
%   f        - 2×1 vector-valued function f(x)
%   J        - 2×2 Jacobian matrix J_f(x)
%
% Dependencies:
%   None (self-contained)
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> [f, J] = f_and_J([1; 3], 1e-6)
%
% ========================================================================

function [f, J] = f_and_J(x, epsilon)

    %% === Extract inputs ===
    x1 = x(1);
    x2 = x(2);
    y = x2 - 3;

    %% === Evaluate nonlinear function f(x) ===
    f1 = (x1 - 1) + y^2;
    f2 = epsilon * y + (3/2) * (x1 - 1) * y + y^2 + y^3;
    f = [f1; f2];

    %% === Compute Jacobian matrix J_f(x) ===
    df1_dx1 = 1;
    df1_dx2 = 2 * y;
    df2_dx1 = (3/2) * y;
    df2_dx2 = epsilon + (3/2) * (x1 - 1) + 2 * y + 3 * y^2;

    J = [df1_dx1, df1_dx2;
         df2_dx1, df2_dx2];
end
