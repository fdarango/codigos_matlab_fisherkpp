function x0 = ic_fkpp(Nm, d)
    Nq = 400;
    x = linspace(0, d, Nq);
    u0 = sin(pi*x/d);

    phi = cell(Nm,1);
    for n = 1:Nm
        phi{n} = sqrt(2/d) * sin(n * pi * x / d);
    end

    x0 = zeros(Nm,1);
    for n = 1:Nm
        numerator = trapz(x, u0 .* phi{n});
        x0(n) = numerator;
end
