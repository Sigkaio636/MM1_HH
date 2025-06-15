function dxdt = HHredu1(t, x)
    % Shifted Nernst equilibrium potentials at mV, supercritical Andronov-Hopf bifurcation
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


    v=x(1); n=x(2); i=x(3);
    dxdt = zeros(3,1);

    % Hypothesis reduction. 
    minf = am(v) ./ (am(v) + bm(v));
    hreg = (0.8882 - 1.041 * n);

    dxdt(1) = (i - gK * n.^4 .* (v - EK) - gNa * minf.^3 .*hreg .* (v - ENa) - gL * (v - EL)) / C ;
    dxdt(2) = an(v).*(1-n) - bn(v).*n;
    dxdt(3) = 0;
end