% ========================================================================
% PHASE 2 — Residual and Jacobian for Nonlinear Robin BVP (FD Method)
% ========================================================================
% File: f_and_J_robin_fd.m
% Author: Noah Pettinato
%
% Summary:
%   Constructs the nonlinear residual F(U) and analytic Jacobian JF(U) for
%   the steady-state heat equation with radiation-type Robin boundary
%   conditions, using a one-sided finite difference (FD) discretization.
%   The Robin flux at x = 1 is enforced via a backward difference, yielding
%   a first-order baseline to compare against the second-order CCFD scheme.
%
% Governing PDE:
%       -k uₓₓ = f(x),     u(0) = 0,
%       -k uₓ(1) = α (u(1)⁴ − u*⁴) + g
%
% Key features:
%   • One-sided finite difference enforcement of the nonlinear Robin BC
%   • Analytic Jacobian assembly for Newton’s method
%   • Baseline (1st order) for accuracy comparison with CCFD (ghost point)
%
% Inputs:
%   U       - N×1 vector of interior unknowns [u₁,…,u_N]^T
%   f_fun   - function handle for f(x), evaluated at interior nodes
%   alpha   - radiation coefficient α
%   g       - boundary source term at x = 1
%   u_star  - reference temperature u*
%   h       - grid spacing
%   k       - thermal conductivity
%
% Outputs:
%   F       - N×1 residual vector
%   J       - N×N sparse Jacobian matrix
%
% Dependencies:
%   None (self-contained). Typically called by a Newton solver (e.g., newton_robin_fd.m).
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> [F, J] = f_and_J_robin_fd(U, f_fun, alpha, g, u_star, h, k)
%
% ========================================================================

function [F, J] = f_and_J_robin_fd(U, f_fun, alpha, g, u_star, h, k)

    N = length(U);
    F = zeros(N, 1);
    J = sparse(N, N);

    x = linspace(0, 1, N+2)';  % Includes boundaries
    fx = f_fun(x(2:N+1));      % f(x_j) for interior points

    %% Left boundary at x = 0 (Dirichlet BC: U_0 = 0)
    F(1) = -k * (0 - 2*U(1) + U(2)) / h^2 - fx(1);
    J(1,1) = 2*k / h^2;
    J(1,2) = -k / h^2;

    %% Interior points (j = 2:N-1)
    for j = 2:N-1
        F(j) = -k * (U(j-1) - 2*U(j) + U(j+1)) / h^2 - fx(j);
        J(j,j-1) = -k / h^2;
        J(j,j)   = 2*k / h^2;
        J(j,j+1) = -k / h^2;
    end

    %% Right boundary at x = 1 (Robin BC with backward difference)
    UN   = U(N);
    UNm1 = U(N-1);
    nonlinear = alpha * (UN^4 - u_star^4) + g;
    F(N) = (k / h) * (UNm1 - UN) - nonlinear;
    J(N,N-1) = k / h;
    J(N,N)   = -k / h - 4 * alpha * UN^3;
end