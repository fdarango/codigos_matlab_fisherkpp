% parametros
D = 0.05;
d = 30;
r = 0.05;
alpha = 0.001;
Tf = 365;

% tiempo de ejecucion semi-implicito
tic;
[u_fd, x_fd, t_fd] = semi_implicit_fkpp(500, Tf, D, d, r, alpha, 0);
t1 = toc;

v = 3:5:64;
errors_gal = zeros(size(v));
ts_gal = zeros(size(v));

for idx = 1:length(v)
    Nm = v(idx);
    
    % condiciones iniciales
    x0_gal = ic_fkpp(Nm, d);

    % --- Construir sistema Galerkin
    [N_fun, phi_fun, xgrid] = galerkin_seno_fkpp(Nm, D, d, r, alpha);

    % --- Resolver ODE Galerkin en los mismos tiempos que FD
    options = odeset('RelTol',1e-4,'AbsTol',1e-4);
    tic;
    [T_gal, X] = ode15s(@(t,X) funode(t, X, N_fun), t_fd, x0_gal, options);
    t_gal = toc;
    % --- Reconstruir solución Galerkin u(x,t) en la malla espacial x_fd y tiempos t_fd
    u_gal = zeros(length(t_fd), length(x_fd));

    for i = 1:length(t_fd)
        for j = 1:Nm
            u_gal(i,:) = u_gal(i,:) + X(i,j)*phi_fun{j}(x_fd);
        end
    end

    % eror relativo
     error_matrix = abs(u_fd' - u_gal);
    err_L2 = norm(error_matrix, 'fro') / norm(u_fd, 'fro');
    ts_gal(idx) = t_gal;
    errors_gal(idx) = err_L2;

end

figure; % Abre una nueva figura
hold on;
scatter(ts_gal, errors_gal, 50, 'filled', 'b'); % '50' es el tamaño del marcador, 'filled' lo rellena, 'b' es azul
 % Mantiene el gráfico actual para añadir más elementos

% 3. Etiquetar cada punto
for i = 1:length(ts_gal)
    offset_x = 0.01; % Pequeño desplazamiento en X
    offset_y = 0.01; % Pequeño desplazamiento en Y

    text(ts_gal(i)+offset_x , errors_gal(i)+offset_y , num2str(v(i)), ...
         'VerticalAlignment', 'bottom', ... % Alinea el texto desde la parte inferior
         'HorizontalAlignment', 'left', ...  % Alinea el texto a la izquierda
         'FontSize', 8, ...                % Tamaño de la fuente
         'Color', 'k');                    % Color del texto (negro)
end

% 4. Ajustes del gráfico (opcional)
title('Tiempos de ejecucion vs errores relativos entre métodos');
xlabel('Tiempos de ejecución(s)');
ylabel('Errores relativos');

x1 = 0;
y1 = 1e-2; % Esto es 0.01
 % Por ejemplo, un valor para t2. Ajusta este valor según tu necesidad.
x2 = t1;
y2 = 0;

% 2. Crear un vector con las coordenadas X y otro con las coordenadas Y
% Esto es lo que 'plot' espera: un vector de X's y un vector de Y's
% Los puntos deben estar en el orden (x_inicio, x_fin) y (y_inicio, y_fin)
x_coords = [x1, x2];
y_coords = [y1, y2];

% 3. Dibujar el segmento de línea punteado
plot(x_coords, y_coords, '--k', 'LineWidth', 1.5);
scatter(x1, y1, 80, 'o', 'filled', 'r'); % Punto de inicio: círculo azul relleno
scatter(x2, y2, 80, 'o', 'filled', 'r'); % Punto de fin: círculo rojo relleno

% 5. Etiquetar el primer punto (0, 10e-3)
% Usaremos el valor de 'y1' como etiqueta, formateado para legibilidad
offset_x_start = -3.5; % Pequeño desplazamiento en X
offset_y_start = 0; % Pequeño desplazamiento en Y para el valor de y1

text(x1 + offset_x_start, y1 + offset_y_start, sprintf('Error de referencia (0, %.0e)', y1), ...
     'VerticalAlignment', 'bottom', ...
     'HorizontalAlignment', 'left', ...
     'FontSize', 9, ...
     'FontWeight', 'bold', ...
     'Color', 'r');

% 6. Etiquetar el segundo punto (t2, 0)
% Usaremos el valor de 't2' como etiqueta
offset_x_end = 2; % Desplazamiento negativo para colocar la etiqueta a la izquierda
offset_y_end = 0.008; % Desplazamiento negativo para colocar la etiqueta debajo

text(x2 + offset_x_end, y2 + offset_y_end, sprintf('Tiempo de ejecución semi-implicito %fs',x2), ...
     'VerticalAlignment', 'top', ... % Alinea el texto desde la parte superior (debajo del punto)
     'HorizontalAlignment', 'right', ... % Alinea el texto a la derecha (para que el final esté en x2)
     'FontSize', 9, ...
     'FontWeight', 'bold', ...
     'Color', 'r');

grid on;
hold off;

% --- Función ODE para Galerkin
function dXdt = funode(~, X, N_fun)
        b_vec = N_fun(X);
        dXdt = b_vec(:);
end
