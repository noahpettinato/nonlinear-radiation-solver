% ========================================================================
% PHASE 0 — Linear Steady-State with Robin Boundary Conditions
% ========================================================================
% File: robin.m
% Author: Noah Pettinato
%
% Summary:
%   Compares first-order one-sided (D_−) vs second-order centered (D_0 with
%   ghost points) finite differences for
%       -u''(x) = f(x),  x ∈ (0,1),
%   with Robin boundary conditions:
%       u'(0) + α₀u(0) = g₀,   -u'(1) = α₁u(1) + g₁.
%
%   Uses u(x) = exp(2x) to manufacture f(x) and consistent boundary data,
%   then verifies observed convergence rates in ||·||_∞ and discrete L² norms.
%
% Key choices:
%   • D_−: one-sided BC discretization (1st order)
%   • D_0: ghost-point BC discretization (2nd order)
%   • Error norms: ∞-norm and √h·||·||₂
%
% Inputs:
%   None (parameters defined within the script)
%
% Outputs:
%   Console output of error norms and convergence rates
%   Log–log plots comparing D_− and D_0 convergence in both norms
%
% Dependencies:
%   None external (self-contained)
%
% Usage:
%   >> robin
%
% ========================================================================

function robin()
    % Robin BC parameters
    alpha0 = 2;
    g0 = 0;
    alpha1 = 2;
    g1 = -4 * exp(2);  % from -u'(1) = alpha1 * u(1) + g1 for u(x) = exp(2x)

    %% Exact solution and RHS
    u_exact = @(x) exp(2 * x);
    f_fun   = @(x) -4 * exp(2 * x);  % from -u'' = f

    % Grid sizes
    hvals = [1e-1, 1e-2, 1e-3, 1e-4];
    E_inf_Dm = zeros(size(hvals));
    E_inf_D0 = zeros(size(hvals));
    E_L2_Dm  = zeros(size(hvals));
    E_L2_D0  = zeros(size(hvals));

    fprintf('%8s %15s %15s %15s %15s\n', ...
        'h', 'E_inf D_-', 'E_2 D_-', 'E_inf D_0', 'E_2 D_0');

    for k = 1:length(hvals)
        h = hvals(k);
        M = round(1 / h);
        x = linspace(0, 1, M+1)';

        %% === D_- Method (first-order accurate) ===
        A1 = sparse(M+1, M+1);
        rhs1 = h^2 * f_fun(x);

        for j = 2:M
            A1(j,j-1) = -1;
            A1(j,j)   = 2;
            A1(j,j+1) = -1;
        end

        % Robin BC at x = 0 using forward difference
        A1(1,1) = -1 - h * alpha0;
        A1(1,2) = 1;
        rhs1(1) = h * g0;

        % Robin BC at x = 1 using backward difference
        A1(M+1, M)   = -1;
        A1(M+1, M+1) = 1 + h * alpha1;
        rhs1(M+1) = -h * g1;

        U1 = A1 \ rhs1;

        %% === D_0 Method (second-order ghost point) ===
        A2 = sparse(M+1, M+1);
        rhs2 = h^2 * f_fun(x);

        for j = 2:M
            A2(j,j-1) = -1;
            A2(j,j)   = 2;
            A2(j,j+1) = -1;
        end

        % Left BC at x = 0 using ghost point
        A2(1,1) = 1 + h * alpha0;
        A2(1,2) = -1;
        rhs2(1) = -h * g0 + (h^2 / 2) * f_fun(x(1));

        % Right BC at x = 1 using ghost point
        A2(M+1, M) = -1;            % Coefficient of u_{M-1}
        A2(M+1, M+1) = 1 + h * alpha1;  % Coefficient of u_M
        rhs2(M+1) = -h * g1 + (h^2 / 2) * f_fun(x(M+1));

        U2 = A2 \ rhs2;

        %% === Compute errors
        utrue = u_exact(x);
        err_inf_Dm = norm(U1 - utrue, inf);
        err_L2_Dm  = sqrt(h) * norm(U1 - utrue, 2);
        err_inf_D0 = norm(U2 - utrue, inf);
        err_L2_D0  = sqrt(h) * norm(U2 - utrue, 2);

        E_inf_Dm(k) = err_inf_Dm;
        E_inf_D0(k) = err_inf_D0;
        E_L2_Dm(k)  = err_L2_Dm;
        E_L2_D0(k)  = err_L2_D0;

        fprintf('%8.1e %15.4e %15.4e %15.4e %15.4e\n', ...
            h, err_inf_Dm, err_L2_Dm, err_inf_D0, err_L2_D0);
    end

    %% === Compute and display convergence rates ===
    fprintf('\nConvergence Rates:\n');
    fprintf('%8s %15s %15s %15s %15s\n', ...
        'h', 'Rate_inf D_-', 'Rate_2 D_-', 'Rate_inf D_0', 'Rate_2 D_0');
    
    for k = 1:length(hvals)-1
        h1 = hvals(k);
        h2 = hvals(k+1);
    
        r_inf_Dm = log(E_inf_Dm(k) / E_inf_Dm(k+1)) / log(h1 / h2);
        r_L2_Dm  = log(E_L2_Dm(k)  / E_L2_Dm(k+1))  / log(h1 / h2);
        r_inf_D0 = log(E_inf_D0(k) / E_inf_D0(k+1)) / log(h1 / h2);
        r_L2_D0  = log(E_L2_D0(k)  / E_L2_D0(k+1))  / log(h1 / h2);
    
        fprintf('%8.1e %15.4f %15.4f %15.4f %15.4f\n', ...
            hvals(k), r_inf_Dm, r_L2_Dm, r_inf_D0, r_L2_D0);
    end


    %% === ∞-norm error plot
    figure;
    loglog(hvals, E_inf_Dm, 'bo-', 'LineWidth', 2, 'DisplayName', '\infty-norm Error D_-');
    hold on;
    loglog(hvals, E_inf_D0, 'rx-', 'LineWidth', 2, 'DisplayName', '\infty-norm Error D_0');
    loglog(hvals, hvals, 'k--', 'DisplayName', 'O(h)');
    loglog(hvals, hvals.^2, 'g--', 'DisplayName', 'O(h^2)');
    xlabel('h'); ylabel('Error (∞-norm)');
    title('Convergence in ∞-norm', 'FontSize', 16);
    legend('Location', 'southeast'); grid on;

    %% === L2 error plot
    figure;
    loglog(hvals, E_L2_Dm, 'bo-', 'LineWidth', 2, 'DisplayName', 'L^2 Error D_-');
    hold on;
    loglog(hvals, E_L2_D0, 'rx-', 'LineWidth', 2, 'DisplayName', 'L^2 Error D_0');
    loglog(hvals, hvals, 'k--', 'DisplayName', 'O(h)');
    loglog(hvals, hvals.^2, 'g--', 'DisplayName', 'O(h^2)');
    xlabel('h'); ylabel('Error (L^2 norm)');
    title('Convergence in L^2 norm', 'FontSize', 16);
    legend('Location', 'southeast'); grid on;
end