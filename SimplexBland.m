function [ xopt,B,message, iter, Zielfktnswert] = SimplexBland( A,b,c,Binit,xB)
%function [ xopt,B ,message, iter] = primalSimplex( A,b,c,Binit,xB )
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
% NAMEN DER PERSONEN AUS DER GRUPPE, DATUM

% Toleranz Definieren!(siehe Blatt)
% Eingabefehler abfangen
% Initialisierung
% Einzelnen Schritte des Algorithmus:
% (1) BTRAN:  
% (2) Pricing:
% (3) FTRAN:
% (4) Ratiotest:
% (5) Update:
    

