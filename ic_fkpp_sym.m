function [a0_sym, b0_sym, c0_sym] = ic_fkpp_sym()
    
    syms x d
    % Funciones base senoidales
    phi1 = sqrt(2/d)*sin(pi*x/d);
    phi2 = sqrt(2/d)*sin(2*pi*x/d);
    phi3 = sqrt(2/d)*sin(3*pi*x/d);

    % Condición inicial 
    u0_sym = sin(pi*x/d);

    % Proyección de u0 sobre cada base
    num1 = int(u0_sym * phi1 , x, 0, d);
    num2 = int(u0_sym * phi2 , x, 0, d);
    num3 = int(u0_sym * phi3 , x, 0, d);

    % Coeficientes iniciales
    a0_sym = simplify(num1);
    b0_sym = simplify(num2);
    c0_sym = simplify(num3);
end
