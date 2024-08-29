clc
clear all 
close all
%% AUFGABE 1
% Konstanten
Bessel= [4156549.426 671255.437 4774279.643;
         4157556.078 659598.467 4775198.257;
         4164901.568 668957.449 4767631.082];
GRS80 = [4157055.229 671337.121 4774764.444;
         4158062.001 659678.732 4775683.150;
         4165408.304 669038.839 4768115.115];
APBessel = [4154589.521 6733343.681 4775681.870];

f = 1/1000;
g = 1000;


% % Berechnungen
% y = [];
% for i = 1:3
%     for k = 1:3
%         y_ = GRS80(i,k);
%         y = [y; y_];
%     end
% end
% 
% A = [];
% for i = 1:3
%     A_ = [1*g 0 0 Bessel(i,1)*f  0              Bessel(i,3)*f   -Bessel(i,2)*f;
%           0 1*g 0 Bessel(i,2)*f -Bessel(i,3)*f  0                Bessel(i,1)*f;
%           0 0 1*g Bessel(i,3)*f  Bessel(i,2)*f -Bessel(i,1)*f    0         ];
%     A = [A;A_];
% end
% 
% x = inv(A'*A) * A'* y;
% x(1:3) = x(1:3) * g
% 
% %% Aufgabe 2
% 
% A = [-1  1  0  0  0;
%       0 -1  1  0  0;
%       0  0 -1  1  0;
%       0  0  0 -1  1;
%       1  0  0  0 -1];
% S_y = eye(5)*(0.003^2) ;
% 
% S_alpha = A * S_y * A';
% R = (eye(5) * 1/0.003) * S_alpha * (eye(5) * 1/0.003); %doppelte von dem was sein sollte why?
% 
% %% Aufgabe 3
% 
% % Konstanten
% y = [20.03; 29.94; 49.97; 19.98; 30.03; 50.20];
% s_y = [0.01; 0.02; 0.03; 0.03; 0.02; 0.015];
% s_apriori = 0.03;
% 
% A = [1 0;
%      0 1;
%      1 1;
%      1 0;
%      0 1;
%      1 1];
% 
% 
% 
% %Berechnungen
% % berechnung von P
% S_y = diag(s_y.^2);
% Q_y = s_apriori^-2 * S_y;
% P = inv(Q_y);
% 
% % berechnung von x_dach
% x_dach1 = inv(A'*P*A)*A'*P*y
% x_dach2 = inv(A'*A)*A'*y

% Berechnungen
%erzeugung von y
y = [];
for i = 1:3
    y_ = GRS80(i,:);
    y = [y; y_'];

end

%erzeugung von A
A = [];
for i = 1:3
    A_ = [1 0 0 Bessel(i,1)*f 0 Bessel(i,3)*f -Bessel(i,2)*f;
          0 1 0 Bessel(i,2)*f -Bessel(i,3)*f 0 Bessel(i,1)*f;
          0 0 1 Bessel(i,3)*f Bessel(i,2)*f -Bessel(i,1)*f 0];
    A = [A; A_];
end

%Ausgleichung
x_dach = (inv(A'*A)*A'*y);
x_dach = [x_dach(1:3); x_dach(4:7)*f];
x_dach = [x_dach(1:4); x_dach(5:7)] % Sollte noch einen multiplikator für bogensecunden bekommen
%Neupunkt berechnung 
Translation = x_dach(1:3); 
Rotation = [x_dach(4) -x_dach(7) x_dach(6);
            x_dach(7) x_dach(4) -x_dach(5);
            -x_dach(6) x_dach(5) x_dach(4)];
NP = Translation + Rotation*APBessel'


%% AUFGABE 2
A = [-1  1  0  0  0;
     0 -1  1  0  0;
     0  0 -1  1  0;
     0  0  0 -1  1;
     1  0  0  0 -1];

Q_y = 0.003^2*eye(5);
Sigma_x = A'*Q_y*A

SE = (1/0.003) * eye(5);
C = SE * Sigma_x * SE

%% Aufgabe 3
%konstanten
y = [20.03; 29.94; 49.97; 19.98; 30.03; 50.20];
sigma_y = [0.01; 0.02; 0.03; 0.03; 0.02; 0.15];
sigma_aprio = 0.03;

% berechnung von P
Sigma_y = diag(sigma_y.^2);
P = inv((1/sigma_aprio^2)*Sigma_y)

%Aufstellung A-Matrix 
A = [1 0;
     0 1;
     1 1;
     1 0;
     0 1;
     1 1];

%Ausgleich
x_dach = inv(A'*P*A)*A'*P*y
x_dach = inv(A'*A)*A'*y




