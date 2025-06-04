function [ap_sym, bp_sym, cp_sym] = galerkin_seno_sym()
    % A es \alpha
    syms x a b c D d r A real

    % Base de funciones seno normalizadas
    phi1 = sqrt(2/d)*sin(pi*x/d);
    phi2 = sqrt(2/d)*sin(2*pi*x/d);
    phi3 = sqrt(2/d)*sin(3*pi*x/d);

    % Solución aproximada u en función de coeficientes
    u = a*phi1 + b*phi2 + c*phi3;

    % Segunda derivada respecto x
    uxx = diff(u, x, 2);

    % Definición del Residuo
    RD = D*uxx + r*u*(1 - u) - A*u;

    % Cálculo de proyecciones
    num1 = int(RD * phi1 , x, 0, d);
    num2 = int(RD * phi2 , x, 0, d);
    num3 = int(RD * phi3 , x, 0, d);

    % Funciones a_p, b_p, c_p para EDO: división por denominador
    ap_sym = simplify(num1);
    bp_sym = simplify(num2);
    cp_sym = simplify(num3);
end
