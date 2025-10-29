% ========================================================================
% PHASE 2 — Newton’s Method (Adaptive Damping) for Nonlinear Robin BVP
% ========================================================================
% File: newton_eq20.m
% Author: Noah Pettinato
%
% Summary:
%   Implements a general damped Newton solver with adaptive step size
%   control for solving F(x) = 0. Used to compute steady-state solutions
%   of the nonlinear Robin boundary value problem defined by either the
%   FD or CCFD discretization. At each iteration, damping (βₖ) is reduced
%   until the residual decreases, ensuring convergence near stiff regimes.
%
% Key features:
%   • Analytic Jacobian supplied via function handle F_and_J(x)
%   • Adaptive damping: halves βₖ until ‖F_new‖ < ‖F_old‖
%   • Monitors ill-conditioning via condest(J)
%   • Terminates on ‖F‖₂ < tol or max_iter reached
%
% Inputs:
%   x0        - initial guess vector
%   F_and_J   - function handle returning [F, J] for the system F(x) = 0
%   tol       - convergence tolerance on residual norm ‖F‖₂
%   max_iter  - maximum number of Newton iterations
%
% Outputs:
%   x         - final approximation to steady-state solution
%   num_iter  - number of iterations performed
%   converged - true if ‖F‖₂ < tol before max_iter
%
% Dependencies:
%   f_and_J_robin_fd.m or f_and_J_robin_ccfd.m  (for residual and Jacobian)
%
% Reproducibility:
%   Deterministic (no randomness)
%
% Usage:
%   >> [x, num_iter, converged] = newton_eq20(x0, F_and_J, 1e-10, 100)
%
% ========================================================================

function [x, num_iter, converged] = newton_eq20(x0, F_and_J, tol, max_iter)

    x = x0;
    converged = false;
    [F, ~] = F_and_J(x);
    res_old = norm(F, 2);   % Initial residual

    for k = 1:max_iter
        [F, JF] = F_and_J(x);
        res = norm(F, 2);   % Current residual
        if res < tol
            converged = true;
            break;
        end
        
        % Warn if Jacobian is ill-conditioned
        if condest(JF) > 1e12
            warning('Jacobian ill-conditioned.');
        end

        fprintf('Iter %d: ||F|| = %.2e\n', k, res);
        dx = JF \ F;    % Newton step
        beta_k = 1;  % Initial damping factor (full step)
        x_new = x - beta_k * dx;

        % Evaluate new residual
        [F_new, ~] = F_and_J(x_new);
        res_new = norm(F_new, 2);

        %% Adaptive damping: backtrack until residual decreases
        while res_new > res && beta_k > 1e-4
            beta_k = beta_k / 2;
            x_new = x - beta_k * dx;
            [F_new, ~] = F_and_J(x_new);
            res_new = norm(F_new, 2);
        end

        % Abort if no damping helps
        if res_new >= res
            warning('Damping failed to reduce residual. Stopping.');
            break;
        end

        x = x_new;
        res_old = res_new;
    end
    num_iter = k;
end