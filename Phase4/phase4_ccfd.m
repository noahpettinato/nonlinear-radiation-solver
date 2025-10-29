% ========================================================================
% PHASE 4 — Fully Implicit Time Integrator with CCFD Ghost-Point Robin BC
% ========================================================================
% File: phase4_ccfd.m
% Author: Noah Pettinato
%
% Summary:
%   Solves c u_t = k uₓₓ on (0,1) with a nonlinear radiation-type Robin
%   boundary at x = 1 using a fully implicit (Backward Euler) update in
%   time and a centered finite-difference scheme with a ghost point (CCFD)
%   in space. Each time step is solved by Newton’s method with an analytic
%   Jacobian from f_and_J_phase4_ccfd.m.
%
% Governing PDE and BCs:
%       c u_t = k uₓₓ,     u(0,t) = 0,
%       -k uₓ(1,t) = α (u(1,t)⁴ − u*⁴) + g
%
% Key features:
%   • Fully implicit (Backward Euler) time discretization
%   • CCFD spatial scheme with ghost point enforcing Robin flux at x = 1
%   • Newton nonlinear solve per time step with analytic Jacobian
%   • Plots final-time solution u(x, T_final)
%
% Inputs:
%   None (parameters are defined within the script)
%
% Outputs:
%   Figure: u(x, T_final) after time marching
%   Console: warnings if Newton fails to converge at any step
%
% Dependencies:
%   newton_eq20.m            % adaptive-damped Newton solver
%   f_and_J_phase4_ccfd.m    % residual and Jacobian for CCFD scheme
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> phase4_ccfd
%
% ========================================================================

function phase4_ccfd()

    % --- Physical and numerical parameters -------------------------
    k = 1; c = 1; alpha = 1; g = 0; u_star = 1;
    h = 1/40; N = round(1/h);
    dt = 0.0001; T_final = 0.1; Nt = round(T_final / dt);   % number of time steps
    x = linspace(0, 1, N+2)';    % full grid including boundaries

    % --- Initial condition: includes ghost point on right ----------
    U = sin(pi * x); % includes all points: U_0 to U_{N+1}
    tol = 1e-8; max_iter = 50; % Relaxed tolerance and more iterations
  
    % --- Time stepping loop ----------------------------------------
    for n = 1:Nt
        U_prev = U;
        FJ = @(U_new) f_and_J_phase4_ccfd(U_new, U_prev, alpha, g, u_star, h, k, c, dt);
        [U, ~, conv] = newton_eq20(U_prev, FJ, tol, max_iter);
        if ~conv
            warning('Newton failed at timestep %d', n);
            disp(U);
            break;
        end
    end

    % --- Plot final result -----------------------------------------
    figure; plot(x, U, 'r.-', 'LineWidth', 1.5);
    xlabel('x'); ylabel('u(x, T)'); title('Phase 4: CCFD Method — Final Time Solution'); grid on;
end