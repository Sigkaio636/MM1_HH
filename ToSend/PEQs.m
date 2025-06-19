close all
I = 0;

V_space = linspace(0, 30, 10000);
n_space = linspace(1e-4, 1 - 1e-4, 1000);

an = @(V) (abs(V-10) < 1e-6) .* 0.1 + (abs(V-10) >= 1e-6) .* 0.01 .* (10 - V) ./ (exp(1 - V / 10) - 1);
bn = @(V) 0.125 * exp(-V / 80);
ninf = @(V) an(V)./(an(V)+bn(V));
% ninf = @(V) 1./(1+exp((-53 -V)/15 ));

% Las isoclinas!

V_ninf = zeros(1000, 1);
for i = 1:1000
    eq = @(V) ninf(V) - n_space(i);
    [res, ~, it] = bisection(-200, 200, 1e-8, 1000, eq);
    % [res, ~, it] = newton(0.5, 1e-4, 1000, eq);
    V_ninf(i) = res(end);
end

V_isoclina = zeros(1000, 1);
for i = 1:1000
    x = @(V) [V, n_space(1001 - i), I];
    eq = @(V) HHredu1_V(x(V));
    if i == 1
        [res, ~, it] = secant(-15, 100, 1e-8, 1000, eq);
    else
        [res, ~, it] = secant(V_isoclina(1002 - i) - 0.5, V_isoclina(1002 - i) + 0.5, 1e-8, 1000, eq);
    end
    V_isoclina(1001 - i) = res(end);
end

% El PEQ
eq = @(V) HHredu1_V([V, ninf(V), I]);

testcurve = [];
for v = V_space
    testcurve = [testcurve eq(v)];
end
figure()
hold on; axis on; 
plot(V_space, testcurve)
plot(V_space, 0*testcurve )
hold off;

[V_PEQ_arr, ~, it] = secant(-10, 0, 1e-8, 1000, eq);
disp(it)
V_PEQ = V_PEQ_arr(end)
V_test = V_PEQ;
HHredu1_V([V_test, ninf(V_test), I])


% La Jacobiana
tol = 1e-4;
jacobiana = zeros(2);
jacobiana(1, 1) = diff(@(x) HHredu1_V(x), [V_PEQ, ninf(V_PEQ), I], [tol, 0, 0]);
jacobiana(1, 2) = diff(@(x) HHredu1_V(x), [V_PEQ, ninf(V_PEQ), I], [0, tol, 0]);
jacobiana(2, 1) = diff(@(x) HHredu1_n(x), [V_PEQ, ninf(V_PEQ), I], [tol, 0, 0]);
jacobiana(2, 2) = diff(@(x) HHredu1_n(x), [V_PEQ, ninf(V_PEQ), I], [0, tol, 0]);
jacobiana
eig(jacobiana)

figure()
grid on; hold on; axis tight
ylim([-150, 150])
title("Isoclinas")
xlabel("Proporción de cosas")
ylabel("Voltaje")
legend()
plot(n_space, V_ninf, "DisplayName", "Isoclina dndt = 0")
plot(n_space, V_isoclina, "DisplayName", "Isoclina dVdt = 0")
plot(ninf(V_PEQ(end)), V_PEQ(end), "-o", "DisplayName", "PEQ")

% Draw a trajectory
Vtest = 24.85;
xr0 = [Vtest; ninf(Vtest); I]; % Initial state <- Play changing the intensity
tspan = [0 50];
[tr, xr] = ode45(@HHredu1, tspan, xr0); 
plot(xr(:,2), xr(:,1), "-", "DisplayName", "Trajectory")
plot(xr(1,2), xr(1,1), "-o", "DisplayName", "Start trajectory")


function [xk, ek, it] = bisection(a, b, tol, itmax, f)
    it = 0;
    tolk = Inf;
    xk = (a + b) ./ 2;

    while it < itmax && tolk >= tol
        if f(a) * f(xk(end)) < 0
            b = xk(end);
        else
            a = xk(end);
        end
        xk = [xk (a + b) ./ 2];
        tolk = abs(xk(end) - xk(end-1));
        it = it + 1;
    end

    ek = xk(1:end - 1) - xk(end);
end

function [xk, ek, it] = newton(x0, tol, itmax, f)
    it = 0;
    tolk = Inf;
    xk = x0;
    while it < itmax && tolk >= tol
        xk = [xk, xk(end) - f(xk(end)) / diff(xk(end), f, 1e-8)];
        tolk = abs(xk(end) - xk(end-1));
        it = it + 1;
    end
    ek = xk(1:end - 1) - xk(end);
end

function [xk, ek, it] = secant(a, b, tol, itmax, f)
    it = 0;
    tolk = Inf;
    xk = [a, b];
    fk_1 = f(a); fk_0 = f(b);

    while it < itmax && tolk >= tol
        xk = [xk, xk(end) - (xk(end) - xk(end-1)) ./ (fk_0 - fk_1) .* fk_0];
        tolk = abs(xk(end) - xk(end-1));
        it = it + 1;
        fk_1 = fk_0;
        fk_0 = f(xk(end));
    end

    xk = xk(2:end);
    ek = xk(1:end - 1) - xk(end);
end

function result = diff (f, x0, h)
    result = (f(x0 + h / 2) - f(x0 - h / 2)) ./ norm(h);
end

function dVdt = HHredu1_V(x)
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

    am = @(V) (abs(V-25) < 1e-6) .* 1 + (abs(V-25) >= 1e-6) .* 0.1 * (25 - V) ./ (exp(2.5 - V / 10) - 1);
    bm = @(V) 4 * exp(-V / 18);

    v = x(1);
    n = x(2);
    i = x(3);

    % Hypothesis reduction. 
    minf = am(v) ./ (am(v) + bm(v));
    hreg = (0.8882 - 1.04 * n);

    dVdt = (i - gK * n.^4 .* (v - EK) - gNa * minf.^3 .*hreg .* (v - ENa) - gL * (v - EL)) / C ;
end

function dndt = HHredu1_n(x)
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

    %an = @(V) 0.01 * vtrap(10 - V, 10);
    an = @(V) (abs(V-10) < 1e-6) .* 0.1 + (abs(V-10) >= 1e-6) .* 0.01 .* (10 - V) ./ (exp(1 - V / 10) - 1);
    bn = @(V) 0.125 * exp(-V / 80);
    %am = @(V) 0.1 * vtrap(25 - V, 10);
    am = @(V) (abs(V-25) < 1e-6) .* 1 + (abs(V-25) >= 1e-6) .* 0.1 * (25 - V) ./ (exp(2.5 - V / 10) - 1);
    bm = @(V) 4 * exp(-V / 18);


    v = x(1);
    n = x(2);
    i = x(3);
    dxdt = zeros(3, 1);

    % Hypothesis reduction. 
    minf = am(v) ./ (am(v) + bm(v));
    hreg = (0.8882 - 1.041 * n);

    dndt = an(v).*(1-n) - bn(v).*n;
end
