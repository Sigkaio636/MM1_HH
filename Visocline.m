close all
an = @(V) (abs(V-10) < 1e-6) .* 0.1 + (abs(V-10) >= 1e-6) .* 0.01 .* (10 - V) ./ (exp(1 - V / 10) - 1);
bn = @(V) 0.125 * exp(-V / 80);
ninf = @(V) an(V)./(an(V)+bn(V));

V_space = linspace(-50, 50, 10000);

% El PEQ
eq = @(V, I) HHredu1_V([V, ninf(V), I]);
figure()
hold on; axis on; 

for i = -20:10:40
    testcurve = [];
    for v = V_space
        testcurve = [testcurve eq(v, i)];
    end
    plot(V_space, testcurve, DisplayName=int2str(i))
end
plot(V_space, 0*testcurve, '--', DisplayName="zero")
legend()
hold off;

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
%% Iterate by V* PEQ potential

Vp_space = linspace(-30, 30, 1000);

% Let's indicate manually all equations 
an = @(V) (abs(V-10) < 1e-6) .* 0.1 + (abs(V-10) >= 1e-6) .* 0.01 .* (10 - V) ./ (exp(1 - V / 10) - 1);
bn = @(V) 0.125 * exp(-V / 80);
ninf = @(V) an(V)./(an(V)+bn(V));
Dan =  @(V) (abs(V-10) < 1e-6) .* 0.01/2 + (abs(V-10) >= 1e-6) .* 0.01 .* exp(V/10).*( 10*exp(V/10) - exp(1)*V ) ./ (10 * (exp(1) - exp(V/10)).^2) ; 
Dbn = @(V) -0.125/80 * exp(-V / 80);

hreg = @(n) 0.8882 - 1.04 * n ;
Dhreg = @(n) -1.04 + 0*n;

am = @(V) (abs(V-25) < 1e-6) .* 1 + (abs(V-25) >= 1e-6) .* 0.1 .* (25 - V) ./ (exp(2.5 - V / 10) - 1);
bm = @(V) 4 * exp(-V / 18);
minf = @(V) am(V)./(am(V)+bm(V));
Dam =  @(V) (abs(V-25) < 1e-6) .* 0.1/2 + (abs(V-25) >= 1e-6) .* 0.1 .* exp(V/10).*( 10*exp(V/10) + exp(2.5)*(15-V) ) ./ ( 10*(exp(2.5) - exp(V/10)).^2 ) ; 
Dbm = @(V) -4/18 * exp(-V / 18);
Dminf = @(V) ( Dam(V).*bm(V) - am(V).*Dbm(V) ) ./ ((am(V)+bm(V)).^2);

% check grafically the formules
%figure
%hold on; axis on; grid on; axis tight;
%plot(Vp_space, Dminf(Vp_space), '-')
%plot(Vp_space, (minf(Vp_space+1e-1) - minf(Vp_space-1e-1))/(2e-1), '-')
%hold off;

%% Calculate the intensity to V* be a critical value

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

igV = @(V) gK .* ninf(V).^4 .* (V-EK) + gNa .* minf(V).^3 .* hreg(ninf(V)) .* (V-ENa) + gL .* (V-EL);
figure
hold on; axis on; grid on; axis tight;
plot(Vp_space, igV(Vp_space), '-')
xlabel("V*")
ylabel("Intensity i for V*")
hold off;

%% Calculate manually the entries of the jacobian at PEQ

VdotV = @(V) ( -gK.*ninf(V).^4 -gNa.* 3.*minf(V).^2.*Dminf(V) .*hreg(ninf(V)).*(V-ENa) -gNa.*minf(V).^3.*hreg(ninf(V)) -gL )/C ;
Vdotn = @(V) ( -gK.*4.*ninf(V).^3.*(V-EK) -gNa.*minf(V).^3.*Dhreg(ninf(V)).*(V-ENa)  )/C ;
ndotV = @(V) ( Dan(V).*bn(V) + an(V).*Dbn(V) )./( an(V)+bn(V) ) ;
ndotn = @(V) bn(V) - an(V) ;

TrJ = @(V) VdotV(V) + ndotn(V) ;
detJ = @(V) VdotV(V).*ndotn(V) - Vdotn(V).*ndotV(V) ;
discM = @(V) TrJ(V).^2 - 4.*detJ(V) ;

figure()
hold on; axis on; grid on;
plot(TrJ(Vp_space), detJ(Vp_space), DisplayName="I evolution")
Det_space = linspace( min(detJ(Vp_space)), max(detJ(Vp_space)), 100 );
Tr_space = linspace( -2*sqrt(max(detJ(Vp_space))), 2*sqrt(max(detJ(Vp_space))), 100 );

plot( 0 .* Tr_space, Tr_space.^2 /4, 'k', DisplayName="Trace =0" )
plot( Tr_space, 0 .* Tr_space.^2 /4, 'k', DisplayName="Determinant =0" )
plot( Tr_space, Tr_space.^2 /4, 'k', DisplayName="Discriminant=0" )
legend(Location='southeast')
hold off;

%% Search for different PEQs evolution

DigV = @(V) (igV(V+1e-6)-igV(V-1e-6))/(2e-6);

[res, ~, it] = bisection(14, 18, 1e-8, 1e3, DigV);
it
maxVp = res(end);
[res, ~, it] = bisection(21, 24, 1e-8, 1e3, DigV);
it
minVp = res(end);

minIbranch = igV(minVp); 
maxIbranch = igV(maxVp);

search_beginBranch_Vp = @(V)igV(V)-minIbranch;
[res, ~, it] = bisection(10.5, 11.5, 1e-8, 1e3, search_beginBranch_Vp);
it
beginBranch_Vp = res(end);

search_endBranch_Vp = @(V)igV(V)-maxIbranch;
[res, ~, it] = bisection(25, 26, 1e-8, 1e3, search_endBranch_Vp);
it
endBranch_Vp = res(end);

lowBranch_Vp_space = linspace(min(Vp_space), maxVp, 10000);
midBranch_Vp_space = linspace(maxVp, minVp, 10000);
higBranch_Vp_space = linspace(minVp, max(Vp_space), 10000);

figure()
hold on; axis on; grid on;
plot(TrJ(lowBranch_Vp_space), detJ(lowBranch_Vp_space), DisplayName="lowBranch PEQ")
plot(TrJ(midBranch_Vp_space), detJ(midBranch_Vp_space), DisplayName="midBranch PEQ")
plot(TrJ(higBranch_Vp_space), detJ(higBranch_Vp_space), DisplayName="higBranch PEQ")
plot(TrJ(beginBranch_Vp), detJ(beginBranch_Vp), 'gx', DisplayName="BeginBranch simple PEQ")
plot(TrJ(minVp), detJ(minVp), 'go', DisplayName="BeginBranch double PEQ")
plot(TrJ(endBranch_Vp), detJ(endBranch_Vp), 'rx', DisplayName="EndBranch simple PEQ")
plot(TrJ(maxVp), detJ(maxVp), 'ro', DisplayName="EndBranch double PEQ")

plot( 0 .* Tr_space, Tr_space.^2 /4, 'k', DisplayName="Trace =0" )
plot( Tr_space, 0 .* Tr_space.^2 /4, 'k', DisplayName="Determinant =0" )
plot( Tr_space, Tr_space.^2 /4, 'k', DisplayName="Discriminant=0" )
legend(Location='southeast')
hold off;


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
