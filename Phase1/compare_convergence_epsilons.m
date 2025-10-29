% ========================================================================
% PHASE 1 — Compare Newton Convergence Domains for Pollock’s Equation
% ========================================================================
% File: compare_convergence_epsilons.m
% Author: Noah Pettinato
%
% Summary:
%   Computes and visualizes Newton’s method convergence domains for several
%   fixed ε values in Pollock’s Equation (Eq. 20). Each figure shows the
%   number of iterations required for convergence across a grid of initial
%   guesses (x₁, x₂), revealing how the domain changes as ε → 0⁺.
%
% Key choices:
%   • ε values: [10⁻³, 10⁻⁴, 10⁻⁵, 10⁻⁶]
%   • Convergence test: ||F||∞ < 1e−10, max_iter = 100
%   • Grid: x₁ ∈ [−2, 2], x₂ ∈ [0, 6], resolution 100 × 100
%   • Visualization: heatmap of iteration counts over (x₁, x₂)
%   • Marks known equilibrium x₊ = (1, 3) for reference
%
% Inputs:
%   None (parameters set in this script)
%
% Outputs:
%   Figures showing convergence domains for each ε in eps_vals
%
% Dependencies:
%   newton_eq20.m  % [x, iters, converged] = newton_eq20(x0, ε, β, tol, max_iter)
%   (pollock_F.m / pollock_J.m if used within newton_eq20)
%
% Reproducibility:
%   Deterministic (no randomness). For consistent color scaling across
%   figures, consider: caxis([1, max_iter + 1]) after imagesc.
%
% Usage:
%   >> compare_convergence_epsilons
%
% ========================================================================

function compare_convergence_epsilons()
% compare_convergence_epsilons   Visualizes Newton's convergence domain for multiple epsilon values.

    %% === Parameter setup ===
    eps_vals = [1e-3, 1e-4, 1e-5, 1e-6];
    beta = 1;
    tol = 1e-10;
    max_iter = 100;

    x1_vals = linspace(-2, 2, 100);
    x2_vals = linspace(0, 6, 100);
    [X1, X2] = meshgrid(x1_vals, x2_vals);

    %% === Loop over different epsilon values ===
    for idx = 1:length(eps_vals)
        epsilon = eps_vals(idx);
        iter_count = zeros(size(X1));

        for i = 1:numel(X1)
            x0 = [X1(i); X2(i)];
            [~, num_iter, converged] = newton_eq20(x0, epsilon, beta, tol, max_iter);
            iter_count(i) = converged * num_iter + (~converged) * (max_iter + 1);
        end

        %% === Plot convergence domain for current epsilon ===
        IterGrid = reshape(iter_count, size(X1));

        figure;
        imagesc(x1_vals, x2_vals, IterGrid);
        set(gca, 'YDir', 'normal');
        colormap(jet);
        colorbar;
        xlabel('Initial x_1', 'FontSize', 14);
        ylabel('Initial x_2', 'FontSize', 14);
        title(['Convergence domain (iterations) for \epsilon = ', num2str(epsilon)], 'FontSize', 12);
        %% === Plot known solution ===
        hold on;
        plot(1, 3, 'ko', 'MarkerSize', 8, 'LineWidth', 2);
    end
end
