function dxdt = HHfull2(t, x)
    % Shifted Nernst equilibrium potentials at m
    EK = -12;
    ENa = 120;
    EL = 10.6;
    % Maximal conductances at mS/cm^2
    gK = 10;
    gNa = 120;
    gL = 0.1;
    % Membrane capacitance and currents
    C = 1;  % µF/cm^2

    % Opening and closing rates
    % Helper function to handle removable singularities in gating variables
    vtrap = @(x, y) (abs(x ./ y) < 1e-6) .* (y .* (1 - x ./ y / 2)) + ...
                    (abs(x ./ y) >= 1e-6) .* (x ./ (exp(x ./ y) - 1)) + ...
                    (abs(x ./ y) < 1e-6 & abs(x) < 1e-6) .* (1e-6);  % Avoid extremely small values
    
    ninf = @(V) 1./(1+exp((-53 -V)/15 ));
    minf = @(V) 1./(1+exp((-40 -V)/15 ));
    hinf = @(V) 1./(1+exp((-60 -V)/(-7) ));
    tn = @(V) 6.1 + 4.7 *exp( -(-79/50-V/50).^2 );
    tm = @(V) 0.04 + 0.46 *exp( -(-38/30-V/30).^2 );
    th = @(V) 1.2 + 7.4 *exp( -(-67/20-V/20).^2 );

    v=x(1); n=x(2); m=x(3); h=x(4); i=x(5);
    dxdt = zeros(5,1);
    dxdt(1) = (i - gK * n.^4 .* (v - EK) - gNa * m.^3 .*h .* (v - ENa) - gL * (v - EL)) / C ;
    dxdt(2) = (ninf(v)-n)./tn(v);
    dxdt(3) = (minf(v)-m)./tm(v);
    dxdt(4) = (hinf(v)-h)./th(v);
    dxdt(5) = 0;
end