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
    %Pricing
    c_N = c(N);
    A_N = A(:, N);
    z_N = c_N - A_N.' * y;
    %Check if optimal
    if all(z_N >= -tol)
        message = 'OPTIMAL VALUE FOUND';
        break
    end
    %Pick random j with z_N_j < 0 (should be at least 1 due to above)
    %{
    valid_j = [];
    s_z_N = size(z_N);
    for j_ind = 1:s_z_N(1)
        if z_N(j_ind) < tol
            valid_j = [valid_j; j_ind];
        end
    end
    s_v_j = size(valid_j);
    j = N(valid_j(randi(s_v_j(1))));
    %}
    %funktioniert nicht so ganz also ist hier die gleiche wie bei Bland
    j_ind = 1;
    while z_N(j_ind) >= -tol
        j_ind = j_ind+1;
    end
    j = N(j_ind);
    %
    %FTRAN A_B*w=A(:j)
    w = linsolve(A_B, A(:,j));
    %Check if unbounded
    if all(w <= tol)
        message = 'PROBLEM IS UNBOUNDED';
        break
    end
    %Ratio Test
    i = 0;
    gamma = Inf;
    A_tmp = A_B\A;
    for k = 1:rank(A)
        if w(k) <= 0
            continue;
        end
        gamma_k = x(B(k)) / w(k);
        if gamma_k - gamma < tol
            gamma = gamma_k;
            i = k;
            continue;
        end
        if abs(gamma_k - gamma) <= tol
            %Lexikografischer Vergleich
            lexikographicVec_k = A_tmp(:,k) / w(k);
            disp(lexikographicVec_k)
            lexikographicVec_i = A_tmp(:,i) / w(i);
            disp(lexikographicVec_i)
            for l = 1:size(lexikographicVec_k)
                if lexikographicVec_k(l) < lexikographicVec_i(l)
                    i = k;
                    break
                end
            end
        end
    end
    %Update
    x(B) = x(B) - gamma * w;
    x(j) = gamma;
    N(j_ind) = B(i);
    N = sort(N);
    B(i) = j;
    B = sort(B);
end

%set final values for return
xopt = x;
xopt(N) = 0;
Zielfktnswert = c.'*x;