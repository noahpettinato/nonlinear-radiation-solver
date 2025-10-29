% ========================================================================
% PHASE 4 — Newton’s Method (Adaptive Damping) for Fully Implicit Step
% ========================================================================
% File: newton_eq20.m
% Author: Noah Pettinato
%
% Summary:
%   Damped Newton solver with adaptive step size for systems F(x) = 0.
%   Used at each fully implicit time step in Phase 4 (FD or CCFD Robin BC).
%   Halves the step size until the residual decreases, improving robustness
%   near stiff, highly nonlinear regimes.
%
% Key features:
%   • Analytic Jacobian supplied via F_and_J handle
%   • Backtracking line search on residual norm ‖F‖₂
%   • Terminates on ‖F‖₂ < tol or if damping fails to reduce ‖F‖₂
%
% Inputs:
%   x0        - initial guess vector
%   F_and_J   - function handle returning [F(x), J(x)]
%   tol       - convergence tolerance on residual norm
%   max_iter  - maximum number of Newton iterations
%
% Outputs:
%   x         - final solution vector
%   num_iter  - number of iterations performed
%   converged - true if ‖F‖₂ < tol before hitting max_iter
%
% Dependencies:
%   Typically paired with: f_and_J_phase4_fd.m or f_and_J_phase4_ccfd.m
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
    res_old = norm(F, 2);

    for k = 1:max_iter
        [F, JF] = F_and_J(x);
        res = norm(F, 2);

        % Check convergence
        if res < tol
            converged = true;
            break;
        end

        fprintf('Iter %d: ||F|| = %.2e\n', k, res);
        dx = JF \ F;
        beta_k = 1;      % Start with full step
        x_new = x - beta_k * dx;
        
        % Evaluate new residual
        [F_new, ~] = F_and_J(x_new);
        res_new = norm(F_new, 2);
       
        % Adaptive damping: reduce beta until residual decreases
        while res_new > res && beta_k > 1e-4
            beta_k = beta_k / 2;
            x_new = x - beta_k * dx;            
            [F_new, ~] = F_and_J(x_new);
            res_new = norm(F_new, 2);
        end
        
        % If no improvement after damping, exit
        if res_new >= res
            warning('Damping failed to reduce residual. Stopping.');
            break;
        end
       
        % Accept damped step
        x = x_new;
        res_old = res_new;
    end
    num_iter = k;
end
