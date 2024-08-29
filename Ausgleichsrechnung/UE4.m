clc
clear all 
close all
%% Aufgabe 2
% Konstanten
%Genaugkeiten
sigma_apriori = 0.0007;
sigma_niv = 0.001;
sigma_trig = 0.01;

% Beobachtungen [m, m]
y_ganz = [0.1255 53;
          0.5326 34;
          0.9583 73;
         -0.2438 26;
         -0.4743 48;
          0.4742 48;
          0.2439 26;
         -0.9581 73;
         -0.5327 34;
         -0.1257 53;
          0.1260 24;
          0.5330 53;
          0.9580 43;
         -0.2440 66;
         -0.4750 58];
y = y_ganz(:,1);


%Erzeugung der gewichtung
Sigma_dH_niv = diag(sigma_niv^2 * (1/1000)*y_ganz(1:10,2));
Sigma_dH_trig = diag(sigma_trig^2 * (1/1000)*y_ganz(11:15,2));
Sigma_y = [Sigma_dH_niv, zeros(10,5);
           zeros(5,10) Sigma_dH_trig];

P = inv((1/sigma_apriori^2)*Sigma_y);

%Definition der A-Matrix
A = [1 0 0 0 0;
    -1 1 0 0 0;
    0 -1 1 0 0;
    0 0 -1 1 0;
    0 0 0 -1 1;
    0 0 0 1 -1;
    0 0 1 -1 0;
    0 1 -1 0 0;
    1 -1 0 0 0;
    -1 0 0 0 0;
    1 0 0 0 0;
    -1 1 0 0 0;
    0 -1 1 0 0;
    0 0 -1 1 0;
    0 0 0 -1 1];

%Berechnung der Ausgleichung
x_dach = inv(A'*P*A)* A'*P*y

%% Aufgabe 3
Sigma_phi = [1 0 0 ;
             0 1 0 ;
             0 0 1 ];












% % A- Matrix
% A = [ 1 0 0 0 0;
%      -1 1 0 0 0;
%      0 -1 1 0 0;
%      0 0 -1 1 0;
%      0 0 0 -1 1;
%      0 0 0 1 -1;
%      0 0 1 -1 0;
%      0 1 -1 0 0;
%      1 -1 0 0 0;
%      -1 0 0 0 0;
%      1 0 0 0 0;
%      -1 1 0 0 0;
%      0 -1 1 0 0;
%      0 0 -1 1 0;
%      0 0 0 -1 1];
% 
% % Berechnungen
% %gewichtung
% s_y_niv = sigma_niv * sqrt(y_ganz(1:10,2));
% s_y_trig = sigma_trig * y_ganz(11:15,2);
% 
% Q_y = (1 / sigma_apriori^2)*[diag(s_y_niv.^2) zeros(10,5); 
%                          zeros(5,10) diag(s_y_trig.^2)];
% 
% P = inv(Q_y);
% 
% %Ausgeichung
% x_dach = inv(A'*P*A) * A'*P*y