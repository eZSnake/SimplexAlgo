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
% Eike Rehwald, Magnus Junker 25.01.2024

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
%Zielfktnswert = 0;
B = Binit;
x = zeros([length(c) 1]);
x(B) = xB;
N = 1:length(c);
N(B) = [];

%check dimensions
A_s = size(A);
b_s = size(b);
c_s = size(c);
if A_s(1) ~= b_s(1) || A_s(2) ~= c_s(1)
    error('Dimensionen sind nicht kompatibel');
end
%check rank
if rank(A) < A_s(1)
    error('A hat keinen vollen Zeilerang');
end
%check starting base
if rank(A(:, B)) < A_s(1) || any(A(:, B)\b < zeros(A_s(1), 1))
    error('B ist keine primal zulässige Startbasis');
end

while true
    iter = iter + 1;
    %BTRAN
    A_B = A(:, B);
    c_B = c(B);
    y = linsolve(A_B.', c_B);
    disp(y);
    %Pricing
    c_N = c(N);
    A_N = A(:, N);
    z_N = c_N - A_N.' * y;
    disp(z_N);
    %Check if optimal
    if all(z_N >= -tol)
        disp('OPTIMAL')
        message = 'True';
        break
    end
%Impl Dantzig rule here (random)
    %Pick first j with z_N_j < 0 (should be at least 1 due to above)
    j_ind = 1;
    while z_N(j_ind) >= -tol
        j_ind = j_ind+1;
    end %need to set j to diff val bc right now it's just 1, 2, 3
    j = N(j_ind);
    disp(j);

    %FTRAN A_B*w=A(:j)
    w = linsolve(A_B, A(:,j));
    disp(w);
    %Check if unbounded
    if all(w <= tol)
        disp('UNBOUNDED')
        message = 'False';
        break
    end
%Impl Dantzig rule here (Lexikographisch)
    %Ratio Test
    gamma_arr = Inf(rank(A), 1);
    k = 1;
    while k <= rank(A)
        %A_tmp = A_B\A;
        %gamma_arr(i) = A_tmp(i) / w(i);
        if x(B(k)) > 0 && w(k) > 0
            %gamma_arr = [gamma_arr; x(b(i)) / w(i)];
            gamma_arr(k) = x(B(k)) / w(k);
        end
        k = k+1;
    end
    disp(gamma_arr)
    [gamma,i] = min(gamma_arr);

    %Update
    x(B) = x(B) - gamma * w;
    x(j) = gamma;
    %tmp_n = N(j);
    N(j_ind) = B(i); %maybe find()?
    B(i) = j;
    %is this already replacing indexes?
    disp(x)
    disp('--------')
end
x_s = size(x);
%set rest of x to 0
xopt = zeros(x_s(1), 1);
xopt(B) = x(B);

Zielfktnswert = c.'*x;
disp(xopt);
disp(Zielfktnswert);