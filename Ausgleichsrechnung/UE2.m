clc
clear all 
close all

%% A1-----------------------------------------------------------------------------------------------------------------------------------------
% IMPORT CONSTANTS
                    % a1, a2, a3, a4
mitarbeitertabelle = [5   6   19  5;  
                      5   8   3   8;
                      0   14  4   10;
                      1   19  4   1;
                      8   10  0   1;
                      6   7   7   6;
                      4   9   17  4];
Arbeitszeit = 360;

% Berechnung

y = ones(7,1) * Arbeitszeit;

x = inv(mitarbeitertabelle'*mitarbeitertabelle) * mitarbeitertabelle' * y


%% A2-----------------------------------------------------------------------------------------------------------------------------------------
% Konstanten
t = -4:1:11;
y = [-26.59;-0.32;17.79;27.04;30.43;27.41;22.44;14.64;7.87;0.25;-4.09;-5.32;-0.14;11.92;32.83;62.92];

% Berechnungen
for i = 1: length(t)
    A(i,:) = [t(i)^0;t(i)^1;t(i)^2;t(i)^3];
end

x = inv(A'*A) * A' * y

k = -4:0.01:11;
for i = 1: length(k)
pol(i) = x(1)*k(i)^0 + x(2)*k(i) + x(3)*k(i)^2 + x(4)*k(i)^3;
end
plot(k, pol )