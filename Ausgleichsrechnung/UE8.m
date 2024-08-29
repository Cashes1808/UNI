clc
clear all 
close all
%% Import der Daten
Data = [322.54, 589.43, 156.87;
        349.43, 547.09, 127.56;
        374.04, 582.39, 170.55;
        409.23, 538.34, 160.80];
Neupunkt = [250, 450];

%% Aufgabe 1

%% Aufgabe 2
%beobachtungen
%   y[m]     x[m]
y = [63.533, 64.643;
     94.234, 64.075;
     112.534, 45.398;
     82.874, 44.012];

e0 = zeros(8, 1);

%treking der inkonsistenzen
e_dachs = e0;
%Ausgleichsloop
while 1
    %Bedingungsfunktionen
        % Erste Bedingung s12 - s43 = 0
        f1_e0_1 = ((y(1,2) + e0(1)) - (y(2,2) + e0(2)))^2;
        f1_e0_2 = ((y(1,1) + e0(3)) - (y(2,1) + e0(4)))^2;
        f1_e0_3 = ((y(3,2) + e0(5)) - (y(4,2) + e0(6)))^2;
        f1_e0_4 = ((y(3,1) + e0(7)) - (y(4,1) + e0(8)))^2;
        f1_e0 = sqrt(f1_e0_1 + f1_e0_2) - sqrt(f1_e0_3 + f1_e0_4);
        
        % zweite Bedingung s14 - s23 = 0
        f2_e0_1 = ((y(1,2) + e0(1)) - (y(4,2) + e0(6)))^2;
        f2_e0_2 = ((y(1,1) + e0(3)) - (y(4,1) + e0(8)))^2;
        f2_e0_3 = ((y(2,2) + e0(2)) - (y(3,2) + e0(5)))^2;
        f2_e0_4 = ((y(2,1) + e0(4)) - (y(3,1) + e0(7)))^2;
        f2_e0 = sqrt(f1_e0_1 + f1_e0_2) - sqrt(f1_e0_3 + f1_e0_4);
    
%     %Aufstellen der B-Matrix 
%         %erste Zeile 
%         B_11 = 2*sqrt(f1_e0_1);
%         B_12 = -2*sqrt(f1_e0_1);
%         B_13 = 2*sqrt(f1_e0_2);
%         B_14 = -2*sqrt(f1_e0_2);
%         B_15 = -2*sqrt(f1_e0_3);
%         B_16 = 2*sqrt(f1_e0_3);
%         B_17 = -2*sqrt(f1_e0_4);
%         B_18 = 2*sqrt(f1_e0_4);
%     
%         B_1 = [B_11 B_12 B_13 B_14 B_15 B_16 B_17 B_18];
%         
%         %zweite Zeile 
%         B_21 = 2*sqrt(f2_e0_1);
%         B_22 = -2*sqrt(f2_e0_2);
%         B_23 = 2*sqrt(f2_e0_2);
%         B_24 = -2*sqrt(f2_e0_2);
%         B_25 = -2*sqrt(f2_e0_3);
%         B_26 = 2*sqrt(f2_e0_3);
%         B_27 = -2*sqrt(f2_e0_4);
%         B_28 = 2*sqrt(f2_e0_4);
%         
%         B_2 = [B_21 B_22 B_23 B_24 B_25 B_26 B_27 B_28];
    
%Aufstellen der B-Matrix 
        %erste Zeile 
        B_11 =  (y(1,2) + e0(1) - y(2,2) - e0(2)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_12 = -(y(1,2) + e0(1) - y(2,2) - e0(2)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_13 =  (y(1,1) + e0(3) - y(2,1) + e0(4)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_14 = -(y(1,1) + e0(3) - y(2,1) + e0(4)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_15 = -(y(3,2) + e0(5) - y(4,2) + e0(6)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
        B_16 =  (y(3,2) + e0(5) - y(4,2) + e0(6)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
        B_17 = -(y(3,1) + e0(7) - y(4,1) + e0(8)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
        B_18 =  (y(3,1) + e0(7) - y(4,1) + e0(8)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
    
        B_1 = [B_11 B_12 B_13 B_14 B_15 B_16 B_17 B_18];
        
        %zweite Zeile 
        B_21 =  (y(1,2) + e0(1) - y(4,2) + e0(6)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_22 = -(y(1,2) + e0(1) - y(4,2) + e0(6)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_23 =  (y(1,1) + e0(3) - y(4,1) + e0(8)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_24 = -(y(1,1) + e0(3) - y(4,1) + e0(8)) * 1/ sqrt(f1_e0_1 + f1_e0_2);
        B_25 = -(y(2,2) + e0(2) - y(3,2) + e0(5)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
        B_26 =  (y(2,2) + e0(2) - y(3,2) + e0(5)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
        B_27 = -(y(2,1) + e0(4) - y(3,1) + e0(7)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
        B_28 =  (y(2,1) + e0(4) - y(3,1) + e0(7)) * 1/ sqrt(f1_e0_3 + f1_e0_4);
        
        B_2 = [B_21 B_22 B_23 B_24 B_25 B_26 B_27 B_28];
    
    
        % B - Matrix
        B = [B_1; B_2];
    
    %Berechnen des wiedersprucs
    w = [f1_e0; f2_e0] - B*e0;
    
    %Ausgleichen der Inkonsistenzen
    e_dach = B'*inv(B*B')*w;
    e_dachs = [e_dachs, e_dach];
    
    %Aufstelen eines abbruchkriteriums
    if norm(e_dach) < 1e-5
        break
    end
    e0 = e_dach;

end

plot(e_dachs(1,:))

%% Testaufgabe 



