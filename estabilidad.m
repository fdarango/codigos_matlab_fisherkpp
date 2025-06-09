% Parámetros fijos
D = 0.05;
d = 30;
r = 0.03;
pi_val = pi;
sqrt2 = sqrt(2);
d32 = d^(3/2);

% Rango de alpha a analizar
alpha_vals = linspace(0.0001, 0.3, 300);

% Resultados
a_eq = zeros(size(alpha_vals));
b_eq = zeros(size(alpha_vals));
c_eq = zeros(size(alpha_vals));
stability = zeros(size(alpha_vals)); % 1 = estable, 0 = inestable
eigvals = NaN(3, length(alpha_vals)); % autovalores del Jacobiano

% Valor inicial de búsqueda
x0 = [4, 0, 0.8];

for i = 1:length(alpha_vals)
    alpha = alpha_vals(i);

    % Definir el sistema como función anónima
    system = @(x) [
        x(1)*(r - alpha) - (D*x(1)*pi_val^3 + (8*sqrt2*x(1)^2*d32*r)/3 + ...
            (32*sqrt2*x(2)^2*d32*r)/15 + (72*sqrt2*x(3)^2*d32*r)/35 - ...
            (16*sqrt2*x(1)*x(3)*d32*r)/15)/(d^2*pi_val);
        
        x(2)*(r - alpha) - x(2)*(420*D*pi_val^3 + 448*sqrt2*x(1)*d32*r + ...
            320*sqrt2*x(3)*d32*r)/(105*d^2*pi_val);
        
        x(3)*(r - alpha) - (9*D*x(3)*pi_val^3 - (8*sqrt2*x(1)^2*d32*r)/15 + ...
            (32*sqrt2*x(2)^2*d32*r)/21 + (8*sqrt2*x(3)^2*d32*r)/9 + ...
            (144*sqrt2*x(1)*x(3)*d32*r)/35)/(d^2*pi_val);
    ];

    % Resolver el sistema con fsolve
    options = optimoptions('fsolve','Display','off');
    [sol, ~, exitflag] = fsolve(system, x0, options);

    if exitflag > 0
        a_eq(i) = sol(1);
        b_eq(i) = sol(2);
        c_eq(i) = sol(3);
        
        % Calcular Jacobiano numéricamente
        J = zeros(3,3);
        delta = 1e-6;
        for j = 1:3
            dx = zeros(3,1); dx(j) = delta;
            J(:,j) = (system(sol + dx) - system(sol - dx)) / (2*delta);
        end
        eigenvals = eig(J);
        eigvals(:, i) = eigenvals;
        
        if all(real(eigenvals) < 0)
            stability(i) = 1; % Estable
        end
        
        % Actualizar x0 para la próxima iteración (continuación)
        x0 = sol;
    else
        a_eq(i) = NaN;
        b_eq(i) = NaN;
        c_eq(i) = NaN;
    end
end

% Dibujar diagrama de bifurcación
figure;
plot(alpha_vals, a_eq, 'b', 'DisplayName', 'a*'); hold on;
plot(alpha_vals, b_eq, 'r', 'DisplayName', 'b*');
plot(alpha_vals, c_eq, 'g', 'DisplayName', 'c*');
plot(alpha_vals(stability==1), a_eq(stability==1), 'bo');
plot(alpha_vals(stability==1), b_eq(stability==1), 'ro');
plot(alpha_vals(stability==1), c_eq(stability==1), 'go');
xlabel('\alpha')
ylabel('Equilibrio')
title('Diagrama de bifurcación respecto a \alpha')
legend show
grid on

% Dibujar parte real de los autovalores
figure;
plot(alpha_vals, real(eigvals(1, :)), 'r', 'DisplayName', 'Re(\lambda_1)'); hold on;
plot(alpha_vals, real(eigvals(2, :)), 'g', 'DisplayName', 'Re(\lambda_2)');
plot(alpha_vals, real(eigvals(3, :)), 'b', 'DisplayName', 'Re(\lambda_3)');
yline(0, '--k', 'DisplayName', 'Re = 0');
xlabel('\alpha');
ylabel('Parte real de los autovalores');
title('Estabilidad: parte real de los autovalores del Jacobiano');
legend show;
grid on

% Dibujar parte imaginaria de los autovalores
figure;
plot(alpha_vals, imag(eigvals(1, :)), 'r', 'DisplayName', 'Im(\lambda_1)'); hold on;
plot(alpha_vals, imag(eigvals(2, :)), 'g', 'DisplayName', 'Im(\lambda_2)');
plot(alpha_vals, imag(eigvals(3, :)), 'b', 'DisplayName', 'Im(\lambda_3)');
yline(0, '--k');
xlabel('\alpha');
ylabel('Parte imaginaria de los autovalores');
title('Parte imaginaria de los autovalores del Jacobiano');
legend show;
grid on

