% ========================================================================
% PHASE 3 — Backward Euler Step with FD Robin Boundary (Backward Diff.)
% ========================================================================
% File: phase3_fd_step.m
% Author: Noah Pettinato
%
% Summary:
%   Advances one timestep for c u_t = k uₓₓ using Backward Euler in time
%   and a one-sided finite difference at x = 1 to enforce the nonlinear
%   radiation-type Robin boundary condition. The interior is implicit in
%   time; the boundary flux uses a backward difference evaluated at the
%   previous step (semi-implicit treatment of the nonlinearity).
%
% Governing PDE:
%       c u_t = k uₓₓ,     u(0,t) = 0,
%       -k uₓ(1,t) = α (u(1,t)⁴ − u*⁴) + g
%
% Key features:
%   • Backward Euler (implicit) time stepping for interior nodes
%   • One-sided FD enforcement of the Robin flux at x = 1
%   • Semi-implicit boundary nonlinearity: α (u⁴ − u*⁴) at previous step
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
%   >> u_new = phase3_fd_step(u_old, alpha, g, u_star, h, dt, c, k)
%
% ========================================================================

function u_new = phase3_fd_step(u_old, alpha, g, u_star, h, dt, c, k)

    N = length(u_old);
    A = sparse(N, N);   % System matrix
    b = zeros(N, 1);    % RHS vector

    % --- Assemble linear system for j = 1 to N-1 (interior) ----------
    main_diag = (1/dt) + (2*k)/(c*h^2);
    off_diag = -k/(c*h^2);

    % Interior rows: j = 1 to N-1
    for j = 1:N-1
        A(j,j) = main_diag;
        if j > 1
            A(j,j-1) = off_diag;
        end
        A(j,j+1) = off_diag;
        b(j) = u_old(j)/dt;
    end

    % --- Apply backward difference for Robin BC at j = N -------------
    nonlinear_bc = alpha * (u_old(N)^4 - u_star^4) + g;
    A(N,N-1) = +k / h;
    A(N,N)   = -k / h;
    b(N)     = nonlinear_bc;

    % --- Solve linear system ----------------------------------------
    u_new = A \ b;
end
