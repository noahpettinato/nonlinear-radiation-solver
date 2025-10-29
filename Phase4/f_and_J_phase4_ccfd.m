% ========================================================================
% PHASE 4 — Residual & Jacobian (Fully Implicit, CCFD Ghost-Point Robin BC)
% ========================================================================
% File: f_and_J_phase4_ccfd.m
% Author: Noah Pettinato
%
% Summary:
%   Builds the nonlinear residual F(U) and analytic Jacobian J(U) for a
%   fully implicit time step of c u_t = k uₓₓ on (0,1) using a centered
%   finite-difference scheme with a ghost point at x = 1 to enforce the
%   radiation-type Robin boundary condition. Left boundary uses a
%   second-order treatment of u(0)=0 via a ghost value.
%
% Governing PDE and BCs:
%       c u_t = k uₓₓ,     u(0,t) = 0,
%       -k uₓ(1,t) = α (u(1,t)⁴ − u*⁴) + g
%
% Key features:
%   • Fully implicit (Backward Euler) residual: (c/Δt)(U − U_prev) − k Δ_h U
%   • Ghost point at x=1 from Robin BC: u_{N+1} = u_{N−1} − (2h/k)·[α(u_N⁴ − u*⁴)+g]
%   • Second-order left boundary via u_{−1} = −U₁
%   • Returns F(U) and sparse J(U) for Newton’s method
%
% Inputs:
%   U       - N×1 current interior solution vector
%   U_prev  - N×1 previous-step interior solution
%   alpha   - radiation coefficient α
%   g       - boundary source term
%   u_star  - reference temperature u*
%   h       - grid spacing
%   k       - thermal conductivity
%   c       - heat capacity
%   dt      - time step size Δt
%
% Outputs:
%   F       - N×1 residual vector
%   J       - N×N sparse Jacobian matrix dF/dU
%
% Dependencies:
%   None (self-contained; typically called inside a Newton time stepper)
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> [F, J] = f_and_J_phase4_ccfd(U, U_prev, alpha, g, u_star, h, k, c, dt)
%
% ========================================================================

function [F, J] = f_and_J_phase4_ccfd(U, U_prev, alpha, g, u_star, h, k, c, dt)

    N = length(U); F = zeros(N,1); J = sparse(N,N);

    % Interior nodes: j = 2 to N-1
    for j = 2:N-1
        F(j) = (c/dt)*(U(j) - U_prev(j)) - k*(U(j-1) - 2*U(j) + U(j+1))/h^2;
        J(j,j-1) = -k/h^2; J(j,j) = c/dt + 2*k/h^2; J(j,j+1) = -k/h^2;
    end
    
    % Left boundary (Dirichlet u=0 approximated via ghost value u_{-1} = -U(1))
    % Equivalent to 2nd-order accurate forward difference at j = 1
    F(1) = (c/dt)*(U(1) - U_prev(1)) - k*(-3*U(1) + U(2))/h^2;
    J(1,1) = c/dt + 3*k/h^2; J(1,2) = -k/h^2;    
    
    % Right boundary (nonlinear Robin via ghost point)
    UN = U(N); UNm1 = U(N-1); Nmod = alpha * (UN^4 - u_star^4) + g; 
    
    % Enforce PDE at j = N using ghost point derived from Robin condition
    % u_{N+1} = u_{N-1} - (2h/k) * Nmod  → used in centered stencil
    % Plug into: u_xx ≈ (u_{N-1} - 2u_N + u_{N+1}) / h^2
    % Residual and Jacobian reflect this substitution
    F(N) = (c/dt)*(UN - U_prev(N)) - k*((UNm1 - UN)/h^2 - Nmod/(k*h));
    J(N,N-1) = -k/h^2; J(N,N) = c/dt + k/h^2 + 4*alpha*UN^3/(k*h);
end
