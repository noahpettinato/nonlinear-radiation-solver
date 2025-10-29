% ========================================================================
% PHASE 3 — Time-Dependent Diffusion with Time-Lagged Robin BC (FD vs CCFD)
% ========================================================================
% File: robin_phase3.m
% Author: Noah Pettinato
%
% Summary:
%   Simulates c u_t = k uₓₓ with a radiation-type Robin boundary enforced
%   using explicit time-lagging at x = 1. The interior uses Backward Euler
%   in time. Compares two spatial boundary discretizations:
%   (1) FD: backward difference at x = 1 (first order)
%   (2) CCFD: ghost-point enforcement at x = 1 (second order)
%   Reports ‖u_FD − u_CCFD‖_∞ at T_final and plots final solutions.
%
% Governing PDE:
%       c u_t = k uₓₓ,     u(0,t) = 0,
%       -k uₓ(1,t) = α (u(1,t−Δt)⁴ − u*⁴) + g
%
% Key features:
%   • Backward Euler time-stepping for the interior
%   • Time-lagged nonlinear Robin flux at x = 1
%   • Comparison of FD (1st order) vs CCFD (2nd order) boundary treatment
%   • Final-time plots and infinity-norm difference
%
% Inputs:
%   None (parameters defined within the script)
%
% Outputs:
%   Console: ‖FD − CCFD‖_∞ at T_final for each α
%   Figures: u(x, T_final) for FD and CCFD across α values
%
% Dependencies:
%   phase3_fd_step.m
%   phase3_ccfd_step.m
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> robin_phase3
%
% ========================================================================

function robin_phase3()

    % --- Parameters ----------------------------
    c = 1; k = 1;   % PDE coefficients
    u_star = 1;     % Target steady state
    g = 0;          % Boundary source
    h = 1/40;        % Spatial grid spacing
    dt = 0.0001;     % Time step size
    %     dt = 0.001;

    T_final = 3;
    %    T_final = 0.5;
    alphas = [0.1, 1.0, 10.0];


    % --- Discretization setup ------------------
    N = round(1/h) - 1;
    x = linspace(0, 1, N+2)';         % includes boundaries
    t_steps = round(T_final / dt);      % Number of time steps

    for alpha = alphas
        fprintf('\n--- Simulating alpha = %.1f ---\n', alpha);

        % Initial condition: u(x,0) = x
        u0 = x(2:end-1);

        % Time stepping
        u_fd   = u0;       % FD scheme solution
        u_ccfd = u0;        % CCFD scheme solution

        % --- Time stepping --------------------- 
        for n = 1:t_steps
            u_fd   = phase3_fd_step(u_fd, alpha, g, u_star, h, dt, c, k);
            u_ccfd = phase3_ccfd_step(u_ccfd, alpha, g, u_star, h, dt, c, k);
        end

        % Compare FD vs. CCFD
        fprintf('alpha = %.1f, ||FD - CCFD||_inf = %.3e\n', alpha, norm(u_fd - u_ccfd, inf));

        % --- Plot final solution --------------
        figure;
        plot(x(2:end-1), u_fd, 'bo-', 'DisplayName', 'FD', 'LineWidth', 1.5);
        hold on;
        plot(x(2:end-1), u_ccfd, 'rx-', 'DisplayName', 'CCFD', 'LineWidth', 1.5);
        title(sprintf('Phase 3: Time-lagged Radiation BC (\\alpha = %.1f)', alpha));
        xlabel('x');
        ylabel('u(x, T_{final})');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end
