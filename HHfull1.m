function dxdt = HHfull1(t, x)
    % Shifted Nernst equilibrium potentials at m
    EK = -12;
    ENa = 120;
    EL = 10.6;
    % Maximal conductances at mS/cm^2
    gK = 36;
    gNa = 120;
    gL = 0.3;
    % Membrane capacitance and currents
    C = 1;  % µF/cm^2

    % Opening and closing rates
    % Helper function to handle removable singularities in gating variables
    vtrap = @(x, y) (abs(x ./ y) < 1e-6) .* (y .* (1 - x ./ y / 2)) + ...
                    (abs(x ./ y) >= 1e-6) .* (x ./ (exp(x ./ y) - 1)) + ...
                    (abs(x ./ y) < 1e-6 & abs(x) < 1e-6) .* (1e-6);  % Avoid extremely small values
    
    %an = @(V) 0.01 * vtrap(10 - V, 10);
    an = @(V) 0.01 * (10 - V) ./ (exp(1 - V / 10) - 1);
    bn = @(V) 0.125 * exp(-V / 80);
    %am = @(V) 0.1 * vtrap(25 - V, 10);
    am = @(V) 0.1 * (25 - V) ./ (exp(2.5 - V / 10) - 1);
    bm = @(V) 4 * exp(-V / 18);
    ah = @(V) 0.07 * exp(-V/20);
    bh = @(V) 1 ./ (exp(3 - V / 10) + 1);

    v=x(1); n=x(2); m=x(3); h=x(4); i=x(5);
    dxdt = zeros(5,1);
    dxdt(1) = (i - gK * n.^4 .* (v - EK) - gNa * m.^3 .*h .* (v - ENa) - gL * (v - EL)) / C ;
    dxdt(2) = an(v).*(1-n) - bn(v).*n;
    dxdt(3) = am(v).*(1-m) - bm(v).*m;
    dxdt(4) = ah(v).*(1-h) - bh(v).*h;
    dxdt(5) = 0;
end