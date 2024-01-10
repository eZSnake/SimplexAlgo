function [xopt,B,message, iter, Zielfktnswert] = SimplexDantzig(A,b,c,Binit,xB)
%function [xopt,B ,message, iter] = primalSimplex(A,b,c,Binit,xB)
%
% Primales Simplexverfahren
%
% Input:  A, b, c   - Daten für LP in primaler Standardform
%                      min c'x s.t. Ax=b, x>=0
%         Binit, xB - primal zulaessige Basis, zugehörige Basislösung
%                     (optional)
% Output: xopt      - optimale Lsg
%         B         - zugehörige Basis
%         message   - Information über Optimallösung oder Unbeschraenktheit
%         iter      - Anzahl der Iterationen
%
% Eike Rehwald, DATUM

% Toleranz Definieren!(siehe Blatt)
% Eingabefehler abfangen
% Initialisierung
% Einzelnen Schritte des Algorithmus:
% (1) BTRAN: y^T*A_B=c_B^T / A_B^T*y=c_B
% (2) Pricing: z_B=c_N-A_N^T*y
% (3) FTRAN: 
% (4) Ratiotest: 
% (5) Update: 

tol = 10^-6;
iter = 0;
message = ""; %true or false
Zielfktnswert = 0;
B = Binit;
x = zeros([length(c) 1]);
x(B) = xB;
N = 1:length(c);
N(B) = [];

%check dimensions
A_s = size(A);
b_s = size(b);
c_s = size(c);
if A_s(1) ~= b_s(1) | A_s(2) ~= c_s(1)
    error('Dimensionen sind nicht kompatibel');
end
%check rank
if rank(A) < A_s(1)
    error('A hat keinen vollen Zeilerang');
end
%check starting base

while iter < 3 %impl actual end case here
    iter = iter + 1;
    %BTRAN
    A_B = A(:, B);
    c_B = c(B);
    y = linsolve(A_B.', c_B);
    %Pricing
    c_N = c(N);
    A_N = A(:, N);
    z_N = c_N - A_N.' * y;
    %Check if optimal
    if all(z_N >= tol)
        disp('ERROR')
        return
    end
    %Pick first j with z_N_j < 0 (should be at least 1 due to above)
    j = 1;
    while z_N(j) >= 0
        j = j+1;
    end
    %FTRAN A_B*w=A(:j)
    w = linsolve(A_B, A(:,j));
    %Check if unbounded
    if all(w <= tol)
        return
    end
    %Ratio Test
    gamma_arr = zeros([rank(A) 1]);
    i = 1;
    while i < rank(A)
        A_tmp = A_B\A;
        gamma_arr(i) = A_tmp(i) / w(i);
        i = i+1;
    end
    [~,gamma] = min(gamma_arr);
    %Update
    x(B) = x(B) - gamma * w;
    x(j) = gamma;
    %set rest of x to 0
    N(j) = B(i); %maybe find()?
    B(i) = j;
    %is this already replacing indexes?
end

xopt = x;
Zielfktnswert = -c.'*x;
disp(xopt);
disp(Zielfktnswert);
