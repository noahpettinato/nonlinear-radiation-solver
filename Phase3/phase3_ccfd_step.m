% ========================================================================
% PHASE 3 — Backward Euler Step with CCFD Ghost-Point Robin BC
% ========================================================================
% File: phase3_ccfd_step.m
% Author: Noah Pettinato
%
% Summary:
%   Advances one timestep for c u_t = k u_xx using Backward Euler in time
%   and a centered control finite difference (CCFD) scheme in space with a
%   ghost point at x = 1 to enforce the nonlinear radiation-type Robin BC.
%   The boundary nonlinearity is evaluated at the previous time level
%   (semi-implicit treatment) for robustness.
%
% Governing PDE:
%       c u_t = k uₓₓ,     u(0,t) = 0,
%       -k uₓ(1,t) = α (u(1,t)⁴ − u*⁴) + g
%
% Key features:
%   • Backward Euler (implicit) time update
%   • CCFD spatial discretization with ghost point at x = 1
%   • Semi-implicit Robin flux: α (u⁴ − u*⁴) evaluated at previous step
%   • Builds and solves a sparse linear system A u_new = b each step
%
% Inputs:
%   u_old  - N×1 solution at previous timestep
%   alpha  - radiation coefficient α
%   g      - boundary source term at x = 1
%   u_star - reference temperature u*
%   h      - spatial grid spacing
%   dt     - time step
%   c, k   - PDE coefficients
%
% Outputs:
%   u_new  - N×1 solution at the next timestep
%
% Dependencies:
%   None (self-contained)
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> u_new = phase3_ccfd_step(u_old, alpha, g, u_star, h, dt, c, k)
%
% ========================================================================

function u_new = phase3_ccfd_step(u_old, alpha, g, u_star, h, dt, c, k)

    N = length(u_old);
    A = sparse(N, N);   % System matrix
    b = zeros(N, 1);    % RHS vector

    % --- Assemble linear system for j = 1 to N-1 (interior) ----------
    main_diag = (1/dt) + (2*k)/(c*h^2);
    off_diag  = -k/(c*h^2);

    % Interior rows: j = 1 to N-1
    for j = 1:N-1
        A(j,j) = main_diag;
        if j > 1
            A(j,j-1) = off_diag;
        end
        A(j,j+1) = off_diag;
        b(j) = u_old(j)/dt;
    end

    % --- Apply ghost point BC at j = N -------------------------------
    UN   = u_old(N);
    nonlin = alpha * (UN^4 - u_star^4) + g;

    % Apply final row from full PDE at j = N
    A(N,N-1) = -2*k / (c*h^2);
    A(N,N)   = (1/dt) + 2*k / (c*h^2);
    b(N)     = u_old(N)/dt - (2/(c*h)) * nonlin;

    % --- Solve linear system ----------------------------------------
    u_new = A \ b;
end
