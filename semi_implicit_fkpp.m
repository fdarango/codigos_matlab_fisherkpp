function [u, x, t] = semi_implicit_fkpp(N, Tf, D, d, r, alpha, g)
    % Tf es el tiempo final
    
    % malla espacial
    delta_x = d / (N + 1);
    
    % Condicion de estabilidad
    if alpha < r
        max_delta_x = sqrt(4 * D / (r - alpha));
        if delta_x > max_delta_x
            N = ceil(d / max_delta_x) - 1;
            delta_x = d / (N + 1);
            fprintf('Adjusted N=%d to satisfy stability (delta_x=%.4f)\n', N, delta_x);
        end
    end
    
    % paso temporal igual a paso espacial
    delta_t = delta_x; 
    Nt = round(Tf / delta_t);

    t = linspace(0, Tf, Nt + 1);
    x = linspace(0, d, N + 2); 
    
    % condicion inicial
    u0 = @(x) sin(pi * x / d);
    u = zeros(N + 2, Nt + 1);
    u(:, 1) = u0(x(:));
    
    % condiciones de frontera
    u(1, :) = 0;
    u(end, :) = 0;
    
    % matriz tridiagonal
    theta = D * delta_t / delta_x^2;
    main_diag = (1 + 2 * theta + alpha * delta_t) * ones(N, 1);
    off_diag = -theta * ones(N - 1, 1);
    A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
    
    % solucion para cada tiempo
    for n = 1:Nt
    % Reacción semi-implícita 
        reaction_coeff = r * delta_t * (1 - u(2:end-1, n));
        A_modified = A - diag(reaction_coeff);
    
    % Vector b: u^n + contribuciones de frontera
        b = u(2:end-1, n);
        b(1) = b(1) + theta * u(1, n+1);     % u(1, n+1) = 0
        b(end) = b(end) + theta * u(end, n+1); % u(end, n+1) = 0
    
    % Resolver el sistema tridiagonal
        u(2:end-1, n+1) = A_modified \ b;
    end
    
    % grafica
    if g == 1
        [X, T] = meshgrid(x, t);
        surf(X, T, u');
        colormap turbo;
        colorbar;
        shading interp;
        xlabel('x');
        ylabel('t');
        zlabel('u(x,t)');
        title(['Semi-Implicit: D=', num2str(D), ', d=', num2str(d), ', r=', num2str(r), ', \alpha=', num2str(alpha)]);
    end
end