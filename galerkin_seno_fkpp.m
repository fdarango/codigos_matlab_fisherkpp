function [N_fun, phi_fun, xgrid] = galerkin_seno_fkpp(Nm, D, d, r, alpha)
    
    % Nm = numero de funciones base

    % definicion variables simbolicas
    syms x
    a = sym('a', [Nm, 1]);
    assume(a, 'real')

    % definicion bases seno
    phi = sym(zeros(Nm,1));
    for n = 1:Nm
        phi(n) = sqrt(2/d) * sin(n*pi*x/d);
    end

    % reconstruccion de la EDP
    u = sum(a .* phi);
    uxx = diff(u, x, 2);
    RD = D*uxx + r*u*(1 - u) - alpha*u;

    % funcion matlab de RD
    RD_fun = matlabFunction(RD, 'Vars', {x, a});

    % funciones matlab de las bases
    phi_fun = cell(Nm,1);
    for n = 1:Nm
        phi_fun{n} = matlabFunction(phi(n), 'Vars', {x});
    end

    % definicion intervalo de integracion
    xgrid = linspace(0, d, 400);

    % Función de proyecciones del residuo
    N_fun = @(a_val) arrayfun(@(i) trapz(xgrid, RD_fun(xgrid, a_val(:)).*phi_fun{i}(xgrid)), 1:Nm);
end