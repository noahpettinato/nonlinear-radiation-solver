% ========================================================================
% PHASE 2 — Convergence Study for Nonlinear Robin BVP (FD vs CCFD)
% ========================================================================
% File: robin_phase2_convergence.m
% Author: Noah Pettinato
%
% Summary:
%   Evaluates and compares the convergence behavior of two steady-state
%   discretizations for the nonlinear Robin boundary value problem:
%   (1) one-sided finite difference (FD, first-order) and
%   (2) centered control finite difference (CCFD, second-order).
%   Errors are measured against a manufactured reference solution for
%   multiple grid spacings h and radiation coefficients α.
%
% Governing PDE:
%       -k uₓₓ = f(x),     u(0) = 0,
%       -k uₓ(1) = α (u(1)⁴ − u*⁴) + g
%
% Key features:
%   • Compares FD vs CCFD accuracy under varying α and h
%   • Reference solution computed from steady linear slope A_ref
%   • Newton’s method (adaptive damping) used for nonlinear solves
%   • Produces log–log error plots verifying first- and second-order rates
%
% Inputs:
%   None (parameters defined within the script)
%
% Outputs:
%   Log–log error plots of ||U − u_ref||∞ vs h for each α
%   Console output of convergence results
%
% Dependencies:
%   newton_eq20.m          % adaptive Newton solver
%   f_and_J_robin_fd.m     % FD residual and Jacobian
%   f_and_J_robin_ccfd.m   % CCFD residual and Jacobian
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> robin_phase2_convergence
%
% ========================================================================

function robin_phase2_convergence()
    f_fun = @(x) 0*x;   % Zero source
    u_star = 1;          % Radiation target
    k = 1;              % Conductivity
    tol = 1e-10;         % Newton tolerance
    max_iter = 100;
    g = 0;              % Source term in BC
    alphas = [0.1, 1.0, 10];  % Test cases
    hs = 1 ./ [10, 20, 40, 80];         % Grid spacings
    
    err_fd = zeros(size(hs));
    err_ccfd = zeros(size(hs));
    
    %% === Loop over alpha values ===
    for a = alphas
        fprintf('\n===== Order of Convergence Test: alpha = %.1f =====\n', a);
        for i = 1:length(hs)
            h = hs(i);
            N = round(1/h) - 1;
            x = linspace(0, 1, N+2)';
            U0 = x(2:end-1); % Initial guess

            % Steady reference slope A solves: -k A = alpha(A^4 - u_star^4)
            A_ref = fzero(@(A) -k*A - a*(A^4 - u_star^4), 1);
            u_ref = A_ref * x(2:end-1);

            % FD
            [U_fd, ~, ~] = newton_eq20(U0, @(U) f_and_J_robin_fd(U, f_fun, a, g, u_star, h, k), tol, max_iter);
            err_fd(i) = norm(U_fd - u_ref, inf);

            % CCFD
            [U_ccfd, ~, ~] = newton_eq20(U0, @(U) f_and_J_robin_ccfd(U, f_fun, a, g, u_star, h, k), tol, max_iter);
            err_ccfd(i) = norm(U_ccfd - u_ref, inf);
        end

        %% === Plot error vs h ===
        figure;
        loglog(hs, err_fd, 'bo-', 'LineWidth', 1.5, 'DisplayName', 'FD Error');
        hold on;
        loglog(hs, err_ccfd, 'rx-', 'LineWidth', 1.5, 'DisplayName', 'CCFD Error');
        
        % Reference lines for order
        loglog(hs, hs * err_fd(1)/hs(1), 'b--', 'DisplayName', 'O(h) Reference');
        loglog(hs, hs.^2 * err_ccfd(1)/hs(1)^2, 'r--', 'DisplayName', 'O(h^2) Reference');

        xlabel('h');
        ylabel('Error ||U - u_{ref}||_{\infty}');
        title(sprintf('Convergence Plot (alpha = %.1f)', a));
        legend('Location', 'northwest');
        grid on;
        hold off;
    end
end
