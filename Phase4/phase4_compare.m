% ========================================================================
% PHASE 4 — Fully Implicit FD vs CCFD (Final-Time Comparison)
% ========================================================================
% File: phase4_compare.m
% Author: Noah Pettinato
%
% Summary:
%   Compares fully implicit time integrations of c u_t = k uₓₓ with a
%   nonlinear radiation-type Robin boundary at x = 1 using two spatial
%   boundary discretizations:
%   (1) FD: backward (one-sided) difference at x = 1
%   (2) CCFD: centered scheme with ghost point at x = 1
%   Both schemes use Newton’s method each time step. Plots final solutions
%   and prints the infinity-norm difference ‖U_FD − U_CCFD‖_∞.
%
% Governing PDE and BCs:
%       c u_t = k uₓₓ,     u(0,t) = 0,
%       -k uₓ(1,t) = α (u(1,t)⁴ − u*⁴) + g
%
% Key features:
%   • Fully implicit (Backward Euler) time discretization
%   • FD vs CCFD boundary enforcement at x = 1
%   • Newton solve per time step with analytic Jacobians
%   • Final-time plot and infinity-norm comparison
%
% Inputs:
%   None (parameters defined within the script)
%
% Outputs:
%   Figure: u(x, T_final) for FD and CCFD
%   Console: ‖FD − CCFD‖_∞ at final time
%
% Dependencies:
%   newton_eq20.m
%   f_and_J_phase4_fd.m
%   f_and_J_phase4_ccfd.m
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> phase4_compare
%
% ========================================================================

function phase4_compare()

    % --- Physical and numerical parameters -------------------------
    k = 1; c = 1; alpha = 10; g = 0; u_star = 1;
    h = 1/40; N = round(1/h) - 1;
    x_full = linspace(0, 1, N + 2)';  % full grid with boundaries
    x = x_full(2:end-1);              % interior points

    dt = 0.05; T_final = 0.1; Nt = round(T_final / dt); % number of time steps
    tol = 1e-8; max_iter = 50;

    % --- Initial condition (interior only) -------------------------
    U_fd   = sin(pi * x);
    U_ccfd = sin(pi * x);

    % --- Time stepping: FD scheme ----------------------------------
    for n = 1:Nt
        U_prev = U_fd;
        FJ = @(U_new) f_and_J_phase4_fd(U_new, U_prev, alpha, g, u_star, h, k, c, dt);
        [U_fd, ~, conv] = newton_eq20(U_prev, FJ, tol, max_iter);
        if ~conv, warning('FD: Newton failed at timestep %d', n); break; end
    end

    % --- Time stepping: CCFD scheme --------------------------------
    for n = 1:Nt
        U_prev = U_ccfd;
        FJ = @(U_new) f_and_J_phase4_ccfd(U_new, U_prev, alpha, g, u_star, h, k, c, dt);
        [U_ccfd, ~, conv] = newton_eq20(U_prev, FJ, tol, max_iter);
        if ~conv, warning('CCFD: Newton failed at timestep %d', n); break; end
    end

    % --- Plot final solutions --------------------------------------
    figure;
    plot(x, U_fd, 'b.-', 'DisplayName', 'FD', 'LineWidth', 1.5); hold on;
    plot(x, U_ccfd, 'r.-', 'DisplayName', 'CCFD', 'LineWidth', 1.5);
    xlabel('x'); ylabel('u(x, T)'); grid on;
    title('Phase 4: FD vs CCFD (Final Time Solution)');
    legend('Location', 'best');

    % --- Compute and print infinity norm error ---------------------
    err_inf = norm(U_fd - U_ccfd, inf);
    fprintf('|| FD - CCFD ||_inf = %.3e\n', err_inf);
end
