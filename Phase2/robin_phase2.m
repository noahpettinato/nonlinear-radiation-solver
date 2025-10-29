% ========================================================================
% PHASE 2 — Nonlinear Robin BVP Solver (FD vs CCFD Comparison)
% ========================================================================
% File: robin_phase2.m
% Author: Noah Pettinato
%
% Summary:
%   Executes Newton’s method for both one-sided finite difference (FD) and
%   centered control finite difference (CCFD) discretizations of the
%   nonlinear steady-state heat equation with radiation-type Robin
%   boundary conditions. Compares convergence behavior and resulting
%   steady-state solutions for varying radiation coefficients α.
%
% Governing PDE:
%       -k uₓₓ = f(x),     u(0) = 0,
%       -k uₓ(1) = α (u(1)⁴ − u*⁴) + g
%
% Key features:
%   • Solves the nonlinear BVP using both FD (1st order) and CCFD (2nd order)
%   • Newton solver with analytic Jacobians from f_and_J_robin_fd / ccfd
%   • Compares iteration counts and convergence status for multiple α
%   • Produces log–log plots of |u(x)| for FD vs CCFD results
%
% Inputs:
%   None (parameters defined within the script)
%
% Outputs:
%   Console summary of convergence success and iteration counts
%   Log–log plots of |u(x)| for each α showing FD vs CCFD behavior
%
% Dependencies:
%   newton_eq20.m          % adaptive Newton solver
%   f_and_J_robin_fd.m     % residual and Jacobian (FD)
%   f_and_J_robin_ccfd.m   % residual and Jacobian (CCFD)
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> robin_phase2
%
% ========================================================================

function robin_phase2()
% robin_phase2   Runs Newton solver for FD and CCFD schemes with nonlinear Robin BC

    % Setup
    f_fun = @(x) 0 * x;        % zero source term
    u_star = 1;                % target steady state
    h = 1/20;                  % Grid spacing
    N = round(1/h) - 1;        % Number of interior points
    x = linspace(0, 1, N+2)';  % Grid: x_0 = 0 to x_{N+1} = 1
    k = 1;                     % Thermal conductivity
    tol = 1e-6;                % Newton tolerance
    max_iter = 100;            % Max Newton iterations

    alphas = [0.1, 1.0, 10.0]; % Test cases

    for alpha = alphas
        g = 0;                 % Boundary term
        U0 = x(2:end-1);       % Linear initial guess from 0 to 1

        % Define function handles
        FJ_fd = @(U) f_and_J_robin_fd(U, f_fun, alpha, g, u_star, h, k); 
        FJ_ccfd = @(U) f_and_J_robin_ccfd(U, f_fun, alpha, g, u_star, h, k);

        % Solve with FD method
        fprintf('\n--- Testing alpha = %.1f with FD ---\n', alpha);
        [U_fd, steps_fd, conv_fd] = newton_eq20(U0, FJ_fd, tol, max_iter);
        fprintf('FD: Converged = %d, Steps = %d\n', conv_fd, steps_fd);

        % Solve with CCFD method
        fprintf('\n--- Testing alpha = %.1f with CCFD ---\n', alpha);
        [U_ccfd, steps_ccfd, conv_ccfd] = newton_eq20(U0, FJ_ccfd, tol, max_iter);
        fprintf('CCFD: Converged = %d, Steps = %d\n', conv_ccfd, steps_ccfd);

        % Plot solutions
        figure;
        loglog(x(2:end-1), abs(U_fd), 'bo-', 'LineWidth', 1.5, 'DisplayName', 'FD');
        hold on;
        loglog(x(2:end-1), abs(U_ccfd), 'rx-', 'LineWidth', 1.5, 'DisplayName', 'CCFD');
        title(sprintf('Log-Log Plot: |u(x)| for \\alpha = %.1f', alpha));
        xlabel('x (log scale)');
        ylabel('|u(x)| (log scale)');
        legend('Location', 'best');
        grid on;
        hold off;
    end
end