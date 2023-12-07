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
% (2) Pricing: 
% (3) FTRAN: 
% (4) Ratiotest: 
% (5) Update: 

tol = 10^-6;
iter = 0;
message = ""; %true or false
Zielfktnswert = 0;
xopt = 0;
B = Binit;
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

while iter < 1
    iter = iter + 1;
    %BTRAN
    A_B = A(:, B);
    c_B = c(B);
    y = linsolve(A_B.', c_B);
    %Pricing
    A_N = A(:, N);
    z_N = c_B - A_N.' * y;
    %Check if optimal
    if all(z_N >= tol)
        return
    end
    %Pick random j with z_N_j < 0
    j = 0;
    while z_N(j) >= 0
        j = j+1;
    end
    %FTRAN A_B*w=A(:j)
    w = linsolve(A_B, A(:,j));
    %Ratio Test
    g_t = [];
    i = 1;
    while i < rank(A)
        [g_t, (inv(A_B) * A)(i) / w(i)]
    end
end