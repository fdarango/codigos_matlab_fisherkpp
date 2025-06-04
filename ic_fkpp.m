function x0 = ic_fkpp(Nm, d)
    % intervalo para integrar numericamente
    Nq = 400;
    x = linspace(0, d, Nq);
    
    % condicion inicial
    u0 = sin(pi * x / d);

    % funciones base
    phi = cell(Nm,1);
    for n = 1:Nm
        phi{n} = sqrt(2/d) * sin(n * pi * x / d);
    end

    % calculo condiciones iniciales
    x0 = zeros(Nm,1);
    for n = 1:Nm
        x0(n) = sum(u0 .* phi{n} .* w) / sum(phi{n}.^2 .* w);
    end
end