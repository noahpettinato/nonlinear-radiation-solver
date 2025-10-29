% ========================================================================
% PHASE 1 — Newton Convergence Domain Visualization (Single ε)
% ========================================================================
% File: convergence_plot.m
% Author: Noah Pettinato
%
% Summary:
%   Evaluates and visualizes the convergence behavior of Newton’s method
%   for Pollock’s Equation (Eq. 20) at a fixed ε = 10⁻⁶. Each grid point
%   (x₁, x₂) represents an initial guess; color intensity indicates the
%   number of iterations required for convergence or the maximum iteration
%   limit if Newton fails to converge.
%
% Key choices:
%   • ε = 10⁻⁶, β = 1
%   • Convergence test: ||F||∞ < 1e−10, max_iter = 100
%   • Grid: x₁ ∈ [−2, 2], x₂ ∈ [0, 6], resolution 100 × 100
%   • Visualization: heatmap of iteration counts across (x₁, x₂)
%   • Marks equilibrium x₊ = (1, 3) with a green triangle
%
% Inputs:
%   None (parameters set in this script)
%
% Outputs:
%   Single figure showing the convergence basin for ε = 10⁻⁶
%
% Dependencies:
%   newton_eq20.m  % [x, iters, converged] = newton_eq20(x0, ε, β, tol, max_iter)
%   (pollock_F.m / pollock_J.m if used within newton_eq20)
%
% Reproducibility:
%   Deterministic (no randomness). To keep colors comparable across runs,
%   optionally set caxis([1, max_iter + 1]) after imagesc.
%
% Usage:
%   >> convergence_plot
%
% ========================================================================

function convergence_plot()
% convergence_plot   Visualizes Newton's method convergence domain.

    %% === Grid and parameter setup ===
    x1_vals = linspace(-2, 2, 100);
    x2_vals = linspace(0, 6, 100);
    [X1, X2] = meshgrid(x1_vals, x2_vals);

    epsilon = 1e-6;
    beta = 1;
    tol = 1e-10;
    max_iter = 100;

    iter_count = zeros(size(X1));

    %% === Evaluate convergence across initial guesses ===
    for i = 1:numel(X1)
        x0 = [X1(i); X2(i)];
        [~, num_iter, converged] = newton_eq20(x0, epsilon, beta, tol, max_iter);
        iter_count(i) = converged * num_iter + (~converged) * (max_iter + 1);
    end

    %% === Plot convergence domain ===
    IterGrid = reshape(iter_count, size(X1));

    figure;
    imagesc(x1_vals, x2_vals, IterGrid);
    set(gca, 'YDir', 'normal');
    colormap(jet);
    colorbar;
    xlabel('Initial x_1', 'FontSize', 14);
    ylabel('Initial x_2', 'FontSize', 14);
    title('Convergence domain for Newton''s method (iterations to converge)', 'FontSize', 12);
    %% === Plot known root ===
    hold on;
    plot(1, 3, 'g^', 'MarkerSize', 8, 'LineWidth', 2);
end
