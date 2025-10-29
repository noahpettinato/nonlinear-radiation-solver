% ========================================================================
% PHASE 4 — Fully Implicit Time Integrator with FD Robin Boundary
% ========================================================================
% File: phase4_fd.m
% Author: Noah Pettinato
%
% Summary:
%   Solves the time-dependent diffusion equation c u_t = k uₓₓ with a
%   nonlinear radiation-type Robin boundary at x = 1 using a fully implicit
%   (Backward Euler) scheme in time and a backward (one-sided) difference
%   at the boundary in space. Each time step is solved by Newton’s method
%   with an analytic Jacobian from f_and_J_phase4_fd.m.
%
% Governing PDE and BCs:
%       c u_t = k uₓₓ,     u(0,t) = 0,
%       -k uₓ(1,t) = α (u(1,t)⁴ − u*⁴) + g
%
% Key features:
%   • Fully implicit (Backward Euler) time-stepping
%   • FD enforcement of nonlinear Robin BC at x = 1
%   • Newton solver with analytic Jacobian per time step
%   • Plots final-time solution u(x, T_final)
%
% Inputs:
%   None (parameters defined within the script)
%
% Outputs:
%   Figure: u(x, T_final) for FD scheme
%   Console: warning message if Newton fails at any step
%
% Dependencies:
%   newton_eq20.m
%   f_and_J_phase4_fd.m
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> phase4_fd
%
% ========================================================================

function phase4_fd()

    % --- Parameters and grid setup ---------------------------------
    k = 1; c = 1; alpha = 1; g = 0; u_star = 1;
    h = 1/40; N = round(1/h) - 1; x = linspace(0, 1, N+2)'; % grid including boundary points
    dt = 0.0001; T_final = 0.1; Nt = round(T_final / dt);        % number of time steps

    % --- Initial condition (interior only) -------------------------
    U = sin(pi * x(2:end-1));   % u(x,0) = sin(pi x)

    % --- Newton solver parameters ----------------------------------
    tol = 1e-8; max_iter = 50;

    % --- Time-stepping loop ----------------------------------------
    for n = 1:Nt
        U_prev = U;

        % Define residual/Jacobian for current timestep
        FJ = @(U_new) f_and_J_phase4_fd(U_new, U_prev, alpha, g, u_star, h, k, c, dt);

        % Newton solve
        [U, ~, conv] = newton_eq20(U_prev, FJ, tol, max_iter);
        if ~conv, warning('Newton failed at timestep %d', n); break; end
    end

    % --- Plot final solution ---------------------------------------
    figure; plot(x(2:end-1), U, 'b.-', 'LineWidth', 1.5);
    xlabel('x'); ylabel('u(x, T)'); title('Phase 4: FD Method — Final Time Solution'); grid on;
end
