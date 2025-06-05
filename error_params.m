% --- Parámetros para exploración
alpha_vals = linspace(0, 0.3, 10);
r_vals = linspace(0.002, 0.04, 10);
% --- Inicializar matriz de errores
error_surface = zeros(length(alpha_vals), length(r_vals));

% --- Parámetros fijos
Nm = 3;
D = 0.05;
d = 30;
Tf = 365;
Nx = 200;

for i = 1:length(alpha_vals)
    for j = 1:length(r_vals)
        alpha = alpha_vals(i);
        r = r_vals(j);

        try
            % --- Solución por diferencias finitas
            [u_fd, x_fd, t_fd] = semi_implicit_fkpp(Nx, Tf, D, d, r, alpha, 0);

            % --- Construcción sistema Galerkin
            [N_fun, phi_fun, ~] = galerkin_seno_fkpp(Nm, D, d, r, alpha);
            x0_gal = ic_fkpp(Nm, d);

            % --- Resolver ODEs de Galerkin
            options = odeset('RelTol',1e-4,'AbsTol',1e-4);
            [~, X] = ode15s(@(t,X) funode(t, X, N_fun), t_fd, x0_gal, options);

            % --- Reconstrucción de la solución de Galerkin
            u_gal = zeros(length(t_fd), length(x_fd));
            for ti = 1:length(t_fd)
                for k = 1:Nm
                    u_gal(ti,:) = u_gal(ti,:) + X(ti,k) * phi_fun{k}(x_fd);
                end
            end

            % --- Cálculo del error ||.||_infinito
            error_matrix = abs(u_fd' - u_gal);
            error_surface(i,j) = max(error_matrix(:));

        catch ME
            warning('Error en alpha = %.4f, r = %.4f: %s', alpha, r, ME.message);
            error_surface(i,j) = NaN;
        end
    end
end

% --- Visualización de la superficie de error
[AlphaGrid, RGrid] = meshgrid(alpha_vals, r_vals);
figure;
surf(AlphaGrid, RGrid, error_surface')
xlabel('\alpha')
ylabel('r')
zlabel('Error ||u_{FD} - u_{Galerkin}||_{\infty}')
title(['Superficie de error entre métodos (Nm = ', num2str(Nm), ')'])
shading interp
colormap turbo
colorbar
view(2)  % 

function dXdt = funode(~, X, N_fun)
        b_vec = N_fun(X);
        dXdt = b_vec(:);
end
