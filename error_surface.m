% parametros
Nm = 60; % numero de modos seno
D = 0.05;
d = 30;
r = 0.05;
alpha = 0.001;
Tf = 365; % tiempo final


[u_fd, x_fd, t_fd] = semi_implicit_fkpp(500, Tf, D, d, r, alpha, 0);
% --- Construir sistema Galerkin
[N_fun, phi_fun, xgrid] = galerkin_seno_fkpp(Nm, D, d, r, alpha);
x0_gal = ic_fkpp(Nm, d);

% --- Resolver ODE Galerkin en los mismos tiempos que FD
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[T_gal, X] = ode15s(@(t,X) funode(t, X, N_fun), t_fd, x0_gal, options);
% --- Reconstruir solución Galerkin u(x,t) en la malla espacial x_fd y tiempos t_fd
u_gal = zeros(length(t_fd), length(x_fd));

for i = 1:length(t_fd)
    for j = 1:Nm
        u_gal(i,:) = u_gal(i,:) + X(i,j)*phi_fun{j}(x_fd);
    end
end
error_matrix = abs(u_fd' - u_gal);
% eror relativo
figure;
surf(x_fd, t_fd, error_matrix)
shading interp
colormap turbo
colorbar
xlabel('x')
ylabel('t')
zlabel('Error absoluto')
title(['Error absoluto |u_{FD} - u_{Galerkin}|: D=', num2str(D), ', d=', num2str(d), ', r=', num2str(r), ', \alpha=', num2str(alpha), ', ',num2str(Nm), ' modos'])
view(2)

function dXdt = funode(~, X, N_fun)
        b_vec = N_fun(X);
        dXdt = b_vec(:);
end