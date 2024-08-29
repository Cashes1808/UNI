clc
clear all
close all
%% Aufgabe 1
% Konstanten 
y = [1.4 1.6 1.5 4.8 1.2 3.7 3.6 3.7]'; %1:4 Känguru A 5:8 Känguru B [m]
sigma_apriori = 1; %[m]

%Berechnung von P
Sigma_y = 0.1^2 * eye(8);
Q_y = (1/sigma_apriori^2)*Sigma_y;
P = inv(Q_y)

%Aufstellung von A
A = [0.5 0.5 0.5 0.5*pi 0 0 0 0;
     0 0 0 0 0.5 0.5*pi 0.5*pi 0.5*pi]';

%Ausgleichung
x_dach = inv(A'*P*A)*A'*P*y

%Berechnung y_dach
y_dach = A*x_dach;

%berechnung e_dach
e_dach = y-y_dach
sigma_postA = (e_dach'*P*e_dach)/(8-2)

Q_x_dach = inv(A'*P*A);
Sigma_x_dach = sigma_postA*Q_x_dach
sigmaAB = sqrt(diag(Sigma_x_dach))