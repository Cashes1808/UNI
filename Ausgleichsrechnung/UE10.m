clc
clear all
close all
%% Aufgabe 2
%import der konstanten Werte
%           dx          dy      dz
data = [  364.492,  181.680,  0.8503;
        -  25.743, -352.340, - 5.7027;
        - 338.778,  170.655,   4.8521;
        - 262.700,   26.220, - 5.9940;
          101.808,  207.920, - 5.1442;
           76.075, -144.438,  10.8464];

%           Y[m]    X[m]    Z[m]
Punkte = [-407.44 79.70, 1244.456;
          -42.93, 261.40, 1245.306;
          -68.67, -90.95, 1239.603;
          -144.75 53.48, 1250.450];

y = [data(:,1); data(:,2); data(:,3)];
x0 = [Punkte(:,1); Punkte(:,2); Punkte(:,3)];
% Die ausgleichung mit datumsdefinition ist ein iteratives Verfahren
% bei dem eine DAtumsmatrix D und eine designmatrix A definiert werden
% müssen.
%zunächst muss der Datumsdefekt definiert werden. der ergibt : d = 3
d = 3;
%anschleißend müssen die faktoren des Datumsdeffekts festgestellt werden.In
%dem Fall sind es die translation in x, y und z Richtung.

% Znächst wird die design Matrix konstruiert. 
a1 = [-1 1 0 0;
      0 -1 1 0;
      1 0 -1 0;
      1 0 0 -1;
      0 1 0 -1;
      0 0 1 -1];
a0 = zeros(6,4);

A = [a1 a0 a0;
     a0 a1 a0;
     a0 a0 a1];

%anschließend wird die Bedingungsmatrix D konstuiert
d1 = ones(1,4);
d0 = zeros(1,4);
D = [d1 d0 d0;
     d0 d1 d0;
     d0 d0 d1];

% Auswahl der TeilspurminimierungSpurminimierung

T_alpha = eye(12);
D_alpha = D*T_alpha;

it = 0;
while 1
%III: f_xo berechnen
f_x0 = A*x0;

% IV: dy berechnen
dy = y - f_x0;
%stem(dy) %quicklook


%V: Die Matrizen A und D müssen zusammengeführt werden die ausgleichung
%erfolgt
N = [A'*A   D_alpha'  ;
       D_alpha zeros(3,3)];
n = [A'*dy; zeros(3,1)];

%VI:  Update der größen
dx_dach = inv(N)*n;
x0 = x0 + dx_dach(1:end-3);

%VII: Stopkriterium
if norm(dx_dach) < 1e-4
    break
end
it = it + 1;
end

x_dach = [x0(1:4) x0(5:8) x0(9:12)]

Qx_dach = inv(N'*N);

% Auswahlmatrix
t1_b = [1 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 0 0 1];
t0_b = zeros(4,4);
T_beta = [t1_b, t0_b, t0_b;
          t0_b, t1_b  t0_b;
          t0_b, t0_b, t1_b];

D_beta = D*T_beta;

%berechnung der S-Matrix 
E = null(A'*A);
S_Beta = eye(12) - E*inv(D_beta*E)*D_beta;

dx_alpha = x_dach - Punkte;
dx_alpha = [dx_alpha(:,1); dx_alpha(:,2); dx_alpha(:,3)];

dx_beta = S_Beta*dx_alpha;
x0 = dx_beta;
x_dach_beta = Punkte + [x0(1:4) x0(5:8) x0(9:12)]


