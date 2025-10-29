% ========================================================================
% PHASE 1 — Animate Convergence Domains for Pollock’s Equation (ε → 0⁺)
% ========================================================================
% File: animate_convergence_epsilons.m
% Author: Noah Pettinato
%
% Summary:
%   Visualizes how the convergence basin of Newton’s method for Pollock’s
%   Equation (Eq. 20) evolves as ε decreases toward zero. For each ε, the
%   method is started from a grid of initial guesses (x₁, x₂) and the
%   number of iterations to convergence is recorded, with failures marked
%   by max_iter + 1.
%
% Key choices:
%   • ε values: logspace(1e−1, 1e−6, 30)
%   • Initial guess grid: x₁ ∈ [−2, 2], x₂ ∈ [0, 6], 100 × 100
%   • Convergence threshold: ||F||∞ < 1e−10, max_iter = 100
%   • Visualization: heatmap of iteration counts over (x₁, x₂), with
%     markers for equilibria x₊ = (1, 3) and x₋ computed from
%     η = 1 − √(1 + 2ε), x₋ = (1 − η², 3 + η)
%
% Inputs:
%   None (parameters set in this script)
%
% Outputs:
%   On-screen animation of convergence basins across ε
%
% Dependencies:
%   newton_eq20.m  % [x, iters, converged] = newton_eq20(x0, ε, β, tol, max_iter)
%   (pollock_F.m / pollock_J.m if used within newton_eq20)
%
% Reproducibility:
%   Deterministic (no randomness). To keep colors consistent across frames,
%   optionally set caxis([1, max_iter + 1]) after imagesc.
%
% Usage:
%   >> animate_convergence_epsilons
%
% ========================================================================

function animate_convergence_epsilons()
% animate_convergence_epsilons
% Animate how Newton's convergence domain changes as epsilon -> 0+

    %% === Parameter and grid setup ===
    eps_vals = logspace(-1, -6, 30);
    beta = 1;
    tol = 1e-10;
    max_iter = 100;

    x1_vals = linspace(-2, 2, 100);
    x2_vals = linspace(0, 6, 100);
    [X1, X2] = meshgrid(x1_vals, x2_vals);

    figure;

    %% === Loop over epsilon values and animate ===
    for idx = 1:length(eps_vals)
        epsilon = eps_vals(idx);
        iter_count = zeros(size(X1));

        for i = 1:numel(X1)
            x0 = [X1(i); X2(i)];
            [~, num_iter, converged] = newton_eq20(x0, epsilon, beta, tol, max_iter);
            iter_count(i) = converged * num_iter + (~converged) * (max_iter + 1);
        end

        %% === Plot convergence domain ===
        IterGrid = reshape(iter_count, size(X1));

        imagesc(x1_vals, x2_vals, IterGrid);
        set(gca, 'YDir', 'normal');
        colormap(jet);
        colorbar;
        xlabel('Initial x_1', 'FontSize', 14);
        ylabel('Initial x_2', 'FontSize', 14);
        title(['Convergence Domain for \epsilon = ', num2str(epsilon)], 'FontSize', 14);
        hold on;
        % Compute and plot both roots
        eta = 1 - sqrt(1 + 2 * epsilon);
        x_minus = [1 - eta^2, 3 + eta];
        
        % Plot both roots
        plot(1, 3, 'g^', 'MarkerSize', 8, 'LineWidth', 2);  % x_+
        plot(x_minus(1), x_minus(2), 'x', 'MarkerSize', 8, 'LineWidth', 2, 'Color', [1, 0.5, 0]);  % x_- (new)

        hold off;

        pause(0.5); 
    end
end
