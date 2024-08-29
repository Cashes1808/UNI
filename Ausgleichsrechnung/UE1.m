clc
clear all 
close all
%% IMPORT CONSTANTS
                    % a1, a2, a3, a4
mitarbeitertabelle = [5   6   19  5;  
                      5   8   3   8;
                      0   14  4   10;
                      1   19  4   1;
                      8   10  0   1;
                      6   7   7   6;
                      4   9   17  4];
Arbeitszeit = 360;

Kombinationen = [1 2 3 4; 
                 2 3 4 5; 
                 3 4 5 6; 
                 4 5 6 7; 
                 1 3 5 7];

%% DATA MANIPULATION
com = 1./kombinator(mitarbeitertabelle, Kombinationen(1,:)) * 360;
for i = 1:4
    for k = 1:4
        dauer(i,k) = mitarbeitertabelle(i,k) * com(i);
    end
end





%% Übung aufgabe 2 mit Fremddaten
 data = readtable("Übungsdaten\Actual_consumption_202301010000_202312312359_Day.csv");
 figure
 plot(data.Date, str2double(data.Total_gridLoad__MWh_CalculatedResolutions))
 
%normierung
 load = data.Total_gridLoad__MWh_CalculatedResolutions;
 load = str2double(load(2:50));
 load = load./sum(load);
 figure(Name= "load normalized")
 plot(load)

 %Hilfsvariablen 
 y = load;
 t = (1:numel(y))';
 
 %I:Scätzwerte
 a = 0;
 b = 0.019;
 c = 0.0050;
 d = 7;
 f = 1/2;
 x0 = [a; b; c; d; f];
 
 i = 0;
 while 1
 % II: Funktionales model
 w_ = 2*pi / numel(y);
 phi_ = pi / numel(y);

 f_x0_1 = x0(1)*t + x0(2);
 f_x0_2 = x0(3)* sin(w_* x0(4)* t + x0(5)*phi_);
 
 f_x0 = f_x0_1 + f_x0_2;

 % III: reduktion
 dy = y- f_x0;

 % IV: A-Matrix 
a1 = t;
a2 = ones(numel(y), 1);
a3 =                  sin(w_* x0(4)* t + x0(5)*phi_);
a4 = x0(3)*w_   *t .* cos(w_* x0(4)* t + x0(5)*phi_);
a5 = x0(3)*phi_*      cos(w_* x0(4)* t + x0(5)*phi_);
A = [a1 a2 a3 a4 a5]; 
%rank(A)
 
%V: Ausgleich
dx_dach = inv(A'*A)*A'*dy

% VI: update
x0 = x0 + dx_dach

% VII: Stopt
if norm(dx_dach) < 1e-3
    break
elseif i == 10
    break
end
i = i + 1;
 end

 y_dach = f_x0 + A*dx_dach;
 e_dach = y-y_dach;

 Qx_dach = inv(A'*A);

 R_x_dach = diag(diag(sqrt((1./Qx_dach))))*Qx_dach*diag(diag(sqrt(1./Qx_dach)))


 Qy_dach = A*Qx_dach*A'
 Qe_dach = eye(numel(y)) - Qy_dach

 figure
 subplot(2,1,1) 
 plot(y)
 hold on 
 plot(y_dach)
 hold off
 legend("y", "y_dach")
 subplot(2,1,2)
 stem(e_dach)
 
 
 
 
 
 
 
 
 
 
 
 
 
 %% Funktions
function [combinationmatrix] = kombinator(matrix2analyze, selctorvector)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[rows cols] = size(matrix2analyze);
len = length(selctorvector);
selector = zeros(rows,cols);

for i = 1:len
    k = selctorvector(i);
    selector(k,:) = selector(k,:) + 1;
end
combinationmatrix = diag(selector' * matrix2analyze)';
end

