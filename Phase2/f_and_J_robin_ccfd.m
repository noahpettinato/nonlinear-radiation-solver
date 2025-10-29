% ========================================================================
% PHASE 2 — Residual and Jacobian for Nonlinear Robin BVP (CCFD Method)
% ========================================================================
% File: f_and_J_robin_ccfd.m
% Author: Noah Pettinato
%
% Summary:
%   Constructs the nonlinear residual F(U) and analytic Jacobian JF(U) for
%   the steady-state heat equation with radiation-type Robin boundary
%   conditions, using a centered control finite difference (CCFD) scheme.
%   A ghost point enforces the nonlinear flux at x = 1 with second-order
%   accuracy for use in a Newton steady-state solve.
%
% Governing PDE:
%       -k uₓₓ = f(x),     u(0) = 0,
%       -k uₓ(1) = α (u(1)⁴ − u*⁴) + g
%
% Key features:
%   • CCFD discretization with ghost point at x = 1 for the Robin BC
%   • Analytic Jacobian assembly for Newton’s method
%   • Nonlinear flux coupling via α (u⁴ − u*⁴)
%   • Returns F(U) and JF(U) for a length-N interior grid
%
% Inputs:
%   U       - N×1 vector of interior unknowns [u₁,…,u_N]^T
%   f_fun   - function handle for f(x)
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
%   None (self-contained). Typically called by a Newton solver (e.g., newton_robin_ccfd.m).
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> [F, J] = f_and_J_robin_ccfd(U, f_fun, alpha, g, u_star, h, k)
%
% ========================================================================

function [F, J] = f_and_J_robin_ccfd(U, f_fun, alpha, g, u_star, h, k)

    N = length(U);
    F = zeros(N,1);
    J = sparse(N,N);

    %% Evaluate f(x) at interior points
    x  = linspace(0,1,N+2).';      
    fx = f_fun( x(2:N+1) );

    %% Left Dirichlet BC (u_0 = 0)
    F(1)    = -k*(0 - 2*U(1) + U(2))/h^2 - fx(1);
    J(1,1)  =  2*k/h^2;
    J(1,2)  = -k/h^2;

    %% Interior points (j = 2:N-1)
    for j = 2:N-1
        F(j)      = -k*(U(j-1) - 2*U(j) + U(j+1))/h^2 - fx(j);
        J(j,j-1)  = -k/h^2;
        J(j,j)    =  2*k/h^2;
        J(j,j+1)  = -k/h^2;
    end

    %% Nonlinear Robin BC at x = 1 using ghost point
    UN   = U(N);
    UNm1 = U(N-1);
    nonlin = alpha*(UN^4 - u_star^4) + g;

    % Ghost value: u_{N+1} = u_{N-1} - (2h/k)*nonlinear
    U_Np1 = UNm1 - (2*h/k)*nonlin;

    % Residual at j = N (includes ghost point)
    F(N) = -k*(UNm1 - 2*UN + U_Np1)/h^2 - fx(N);

    % Jacobian contributions at j = N
    J(N,N-1) = -2*k/h^2;                
    J(N,N)   =  2*k/h^2 + 8*alpha*UN^3/h; 
end
