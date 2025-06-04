function [uxt] = solve_galerkin_seno_fkpp(Nm, Tf, D, d, r, alpha, gcoef, gu, tol)

    % para obtener graficas, gcoef = 1, gu = 1
    % Nm: número de funciones base

    % Condiciones iniciales
    x0 = ic_fkpp(Nm, d);

    % Generar funciones de Galerkin
    [N_fun, phi_fun, xgrid] = galerkin_seno_fkpp(Nm, D, d, r, alpha);

    % Resolver sistema ODE
    options = odeset('RelTol', tol, 'AbsTol', tol);
    [T, X] = ode23s(@(t, X) funode(t, X, N_fun), [0 Tf], x0, options);

    % Reconstrucción solución u(x,t)
    Nt = length(T);
    Nx = length(xgrid);
    uxt = zeros(Nt, Nx);
    for i = 1:Nt
        for j = 1:Nm
            uxt(i,:) = uxt(i,:) + X(i,j)*phi_fun{j}(xgrid);
        end
    end

    % Gráfica de los coeficientes
    if gcoef == 1
        figure
        plot(T, X, 'LineWidth', 1.5)
        legend(arrayfun(@(n) sprintf('a_{%d}(t)',n), 1:Nm, 'UniformOutput', false))
        xlabel('t'), ylabel('Coeficientes'), grid on
        title('Coeficientes en el tiempo')
    end

    if gu == 1
        % Solución u(x,t)
        figure
        surf(xgrid, T, uxt)
        shading interp
        colormap turbo
        colorbar
        xlabel('x'), ylabel('t'), zlabel('u(x,t)')
        title('Solución aproximada u(x,t)')
        view(-10,50)
    end

    function dXdt = funode(~, X, N_fun)
        b_vec = N_fun(X);
        dXdt = b_vec(:);
    end

end
