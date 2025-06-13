clear all
close all
format long

%% Exercici pràctic del parcial de 2023, per provar les funcions

% P en funció de e
P = @ (e) 2 * sqrt(1 - e.^2) ./ (e.^3) .* (3 - 2*e.^2) .* asin(e) - 6 ./ (e.^2) .* (1 - e.^2);

% Representació gràfica de P
z = linspace(1e-10, 1 - 1e-10, 2e3 + 1);

figure(1)
plot(z, P(z))
grid on

% 1º apartat: per P = 0.2, trobar e
equation = @(e) P(e) - 0.2;

figure(2)
plot(z, equation(z)) % Dos arrels: entre 0.996 i 0.998; entre 0.596 i 0.6
grid on

disp("BISECTION METHOD")
[xk1, ek1, it1] = bisection(0.596, 0.6, 1e-10, 1000, equation);
xk1(end)
it1
[xk2, ek2, it2] = bisection(0.996, 0.998, 1e-10, 1000, equation);
xk2(end)
it2

disp("CHORD METHOD")
[xk1, ek1, it1] = chord(0.596, 0.6, 1e-10, 1000, equation);
xk1(end)
it1
[xk2, ek2, it2] = chord(0.996, 0.998, 1e-10, 1000, equation);
xk2(end)
it2

disp("SECANT METHOD")
[xk1, ek1, it1] = secant(0.596, 0.6, 1e-10, 1000, equation);
xk1(end)
it1
[xk2, ek2, it2] = secant(0.996, 0.998, 1e-10, 1000, equation);
xk2(end)
it2

disp("NEWTON METHOD")
[xk1, ek1, it1] = newton(0.596, 1e-10, 1000, equation, '0');
xk1(end)
it1
[xk2, ek2, it2] = newton(0.996, 1e-10, 1000, equation, '0');
xk2(end)
it2

disp(equation(xk1(end)))
disp(equation(xk2(end)))


% 2º apartat: trobar el P a partir del qual no hi ha solució.
% La solució es trobar la posició del màxim de la funció i el valor de la
% funció en aquell punt.

h = 1e-8;
derivada = @(e) (P(e + h) - P(e)) / h;

figure(3)
plot(z, derivada(z))
grid on

[xk, ek, it] = newton(0.9, 1e-10, 1000, derivada, '0');
disp(xk(end))
disp(derivada(xk(end)))
disp(P(xk(end)))



%% Les funcions!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Duu a terme el mètode de la bisecció sobre f en l'interval [a, b]
% amb tolerància tol i màxim d'iteracions itmax.
% Retorna:
% - xk: aproximació de l'arrel a cada iteració.
% - ek: error sobre l'arrel trobada a cada iteració.
% - it: nº iteracions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Duu a terme el mètode de Newton-Rhapson sobre f en x0
% amb tolerància tol i màxim d'iteracions itmax.
% Retorna:
% - xk: aproximació de l'arrel a cada iteració.
% - ek: error sobre l'arrel trobada a cada iteració.
% - it: nº iteracions
%
% Opcionalment, es pot inserir df com a derivada de la funció.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xk, ek, it] = newton(x0, tol, itmax, f, df)
    it = 0;
    tolk = Inf;
    xk = x0;

    while it < itmax && tolk >= tol
        if df == '0'
            xk = [xk, xk(end) - f(xk(end)) / diff(xk(end), f, '0')];
        else
            xk = [xk, xk(end) - f(xk(end)) / df(xk(end))];
        end
        tolk = abs(xk(end) - xk(end-1));
        it = it + 1;
    end

    ek = xk(1:end - 1) - xk(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Duu a terme el mètode de la corda sobre f en l'interval [a, b]
% amb tolerància tol i màxim d'iteracions itmax.
% Retorna:
% - xk: aproximació de l'arrel a cada iteració.
% - ek: error sobre l'arrel trobada a cada iteració.
% - it: nº iteracions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xk, ek, it] = chord(a, b, tol, itmax, f)
    it = 0;
    tolk = Inf;
    xk = (a + b) ./ 2;
    fa = f(a); fb = f(b);

    while it < itmax && tolk >= tol
        xk = [xk, xk(end) - (b - a) ./ (fb - fa) .* f(xk(end))];
        tolk = abs(xk(end) - xk(end-1));
        it = it + 1;
    end

    ek = xk(1:end - 1) - xk(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Duu a terme el mètode de la secant sobre f en l'interval [a, b]
% amb tolerància tol i màxim d'iteracions itmax.
% Retorna:
% - xk: aproximació de l'arrel a cada iteració.
% - ek: error sobre l'arrel trobada a cada iteració.
% - it: nº iteracions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Retorna la derivada aproximada de f en x0.
% Utilitza l'aproximació de h donada.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = diff (x0, f, h)
    if h == '0'
        h = 1e-8;
    end
    result = (f(x0 + h) - f(x0)) ./ h;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Ajusta una recta y = a x + b als vectors x, y.
%
%  Torna el coeficient de regressió r, el pendent a i
%  l''ordenada a l''origen b.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,a,b] = reg_lin(x,y)

    if iscolumn(x)
        x = x';
    end
    
    if iscolumn(y)
        y = y';
    end
    
    if size(x,2) ~= size(y,2)
        fprintf('La mida dels dos vectors no es la mateixa al cridar a reg_lin. \n')
        return
    end
    
    n = size(x,2);
    
    ax = sum(x)/n; ay = sum(y)/n; 
    
    ax2 = sum(x.^2)/n; ay2 = sum(y.^2)/n; axy = sum(x.*y)/n;
    
    a = (axy-ax*ay)/(ax2-ax^2);
    c = (axy-ax*ay)/(ay2-ay^2);
    r = sqrt(a*c);
    b = ay-a*ax;
    d = ax-c*ay;
    
    %n,r,a,b,ax,ay,ax2,ay2,axy;

end