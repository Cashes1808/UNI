clc
clear all
close all
format short


%% Festweret und daten
data = readtable('data.txt');
sigma_NH = 0.001; sigma_EH = 0.005;

%% Aufgabenteil a
%Berechnung der höhenanomalie
Zeta = data.EH(1:20)-data.NH(1:20);
sigma_zeta = sqrt(sigma_EH^2+sigma_NH^2);



%darstellung der höhenanomalie

figure;
scatter3(data.X(1:20),data.Y(1:20),data.NH(1:20),40,Zeta,'filled');
text(data.X(1:20)+100,data.Y(1:20)+100,data.NH(1:20)+20,num2str(data.PN(1:20)));
xlabel('x [m]');
ylabel('y [m]');
zlabel('H_N [m]');
c = colorbar;
c.Label.String = 'Höhenanomalie \zeta_i';





Result_A = [data(1:20,1) , table(Zeta)];
writetable(Result_A,'Result_A','FileType','spreadsheet','WriteRowNames',true)

%% Aufgabenteil b
variablen = ['a0'; 'a1'; 'a2'; 'a3'; 'a4'; 'a5'];

%Schwerpunkt der Koordinaten berechnen
x_schwerP = sum(data.X)/length(data.X);
y_schwerP = sum(data.Y)/length(data.Y);

%Reduzierte Koordinaten
x_red = data.X - x_schwerP;
y_red = data.Y - y_schwerP;

% Ausgleichung
% x_dach strucktur [a_0; a_1; a_2; a_3; a_4; a_5]
A = [ones(20,1), y_red(1:20), x_red(1:20), y_red(1:20).*x_red(1:20), y_red(1:20).^2, x_red(1:20).^2];
x_dach = inv(A'*A)*(A'*Zeta);


%berechnung der zusatzangaben
%y_dach(zeta dach)
zeta_dach = A*x_dach;
%e_dach
e_dach = Zeta - zeta_dach;
test = A'*e_dach;
%Berechnung der genauigkeit
sigma_hdach = e_dach'*e_dach/(length(Zeta)-length(x_dach));

Sigma_zeta = sigma_zeta^2*eye(20);
Sigma_adach = sigma_hdach*inv(A'*A);
Sigma_zetadach = A*Sigma_adach*A';
Sigma_edach = Sigma_zeta - Sigma_zetadach;

sigma_adach = sqrt(diag(Sigma_adach));
sigma_zetadach = sqrt(diag(Sigma_zetadach));
sigma_edach = sqrt(diag(Sigma_edach));

Result_B1 = table(variablen, x_dach, sigma_adach);
writetable(Result_B1,'Result_Ausgleichung','FileType','spreadsheet','WriteRowNames',true)

Result_B2 = [data(1:20,1),table(zeta_dach,sigma_zetadach,e_dach,sigma_edach)];
writetable(Result_B2,'Result_Sigmas','FileType','spreadsheet','WriteRowNames',true)

%% Augabenteil C

T = [(1:length(x_dach))', abs(x_dach./sigma_adach)];
alpha = 0.05;
T_quant = tinv(1-alpha/(2*length(T)),length(data.NH(1:20))-length(T));
T_quant_tbl = [T_quant; nan(length(T)-1,1)];
i = 1;



while 1
    %Eliminierungsprozess
    A= A(:,T(:,2)~=min(T(:,2)));
    P= eye(20);
    
    x_dach = (A'*A)\(A'*Zeta);
    x_dach_tbl = [x_dach; nan(length(T)-length(x_dach),1)];

    fprintf('a_%d eliminiert \n',T(T(:,2)==min(T(:,2)),1)-1)
    
    Result_C = table(T(:,2),T_quant_tbl,x_dach_tbl); 
    name = strcat('Result_C',num2str(i));
    writetable(Result_C,name,'FileType','spreadsheet','WriteRowNames',true)
    

    zeta_dach = A*x_dach;
    e_dach = Zeta - zeta_dach;
    
    sigma_hdach = e_dach'*e_dach/(length(Zeta)-length(x_dach));
    
    Sigma_zeta = sigma_zeta^2*eye(20);
    Sigma_adach = sigma_hdach*inv(A'*A);
    Sigma_zetadach = A*Sigma_adach*A';
    Sigma_edach = Sigma_zeta - Sigma_zetadach;
    
    sigma_adach = sqrt(diag(Sigma_adach));
    sigma_zetadach = sqrt(diag(Sigma_zetadach));
    sigma_edach = sqrt(diag(Sigma_edach));
    
    T = [T(T(:,2)~=min(T(:,2)),1) abs(x_dach./sigma_adach)];
    T_quant = tinv(1-alpha/(2*length(T)),length(data.NH)-length(T));
    T_quant_tbl = [T_quant; nan(length(T)-1,1)];
    
    i = i+1;
    if min(T(:,2))<=T_quant
        break
    end
   
end

Result_C = table(x_dach_tbl); 
name = strcat('Result_C',num2str(i));
writetable(Result_C,name,'FileType','spreadsheet','WriteRowNames',true)

%% Aufgabenteil D

% zeta =a(1) + a(2)*y + a(3)*x + a(4)*y*x + a(5)*y^2 +a(6)*x^2;

a = zeros(6,1);
a(T(:,1)) = x_dach;

x_neu = x_red(20+1:end);
y_neu = y_red(20+1:end);
h_neu = data.EH(20+1:end);

A_neu = [ones(length(x_neu),1), y_neu, x_neu, y_neu.*x_neu, y_neu.^2, x_neu.^2];
zeta_neu = A_neu*a;
H_Neu = h_neu-zeta_neu;

Result_D = table(zeta_neu,H_Neu);
writetable(Result_D,'Result_D','FileType','spreadsheet','WriteRowNames',true)


% Ortsbezogene Darstellung der Höhenanomalie
figure;
scatter3(data.X,data.Y,[data.EH(1:20);H_Neu],40,[Zeta;zeta_neu],'filled');
text(data.X+100,data.Y+100,[data.NH(1:20);H_Neu]+20,num2str(data.PN));
xlabel('x [m]');
ylabel('y [m]');
zlabel('H_N [m]');
cb = colorbar;
cb.Label.String = 'Höhenanomalie \zeta_i';


%% Aufgabenteil E

% H_Neu = h_neu - zeta_neu = h_neu a(1) + a(2)*y + a(3)*x + a(4)*y*x + a(5)*y^2 +a(6)*x^2;
% sigma_HNeu = sqrt(sigma_h^2 + sigma_zetaneu^2)
% F_zeta = [dzeta/da_i] = A_neu
Sigma_a = zeros(6,6);
Sigma_a(T(:,1),T(:,1)) = Sigma_adach;
sigma_zetaneu = sqrt(diag(A_neu*Sigma_a*A_neu'));

sigma_HNeu = sqrt(sigma_EH.^2 + sigma_zetaneu.^2);

Result_E = table(sigma_zetaneu,sigma_HNeu);
writetable(Result_E,'Result_E','FileType','spreadsheet','WriteRowNames',true)
