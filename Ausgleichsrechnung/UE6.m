clc
clear all 
close all
%% Aufgabe 1
% Konstanten
%Pos: Länge(lam)[°] Breite(phi)[°] Höhe(h)[°]
Pos = [9.19369 48.80613 262;
       9.19784 48.80574 259;
       9.20012 48.80407 256;
       9.20301 48.80098 247;
       9.20097 48.80136 251;
       9.19808 48.80275 253;
       9.19323 48.80404 254];

%zu interpolierenden Punkte
%P_int: Länge(lam)[°] Breite(phi)[°]
P_int = [9.19905 48.80335;
         9.19414 48.80504];

%Definition der Beobachtungen y
y = Pos(:,3);

% Konstruktion der A - Matrix 
for i = 1:size(Pos,1)
    for k = 1:1:size(Pos,1)
        delta_lambda= (Pos(i,1) - Pos(k,1))^2;
        delta_phi = (Pos(i,2) - Pos(k,2))^2;

        r = sqrt(delta_lambda + delta_phi);
        %festlegung für r = 0 
        if r == 0
            PHI(i,k) = 0;
        else
            PHI(i,k) = r^2 * log(r);
        end
    end
end

A = [ones(size(Pos,1),1) Pos(:,1) Pos(:,2) PHI];

% Konstruktion der bedingungsmatrizen
D_ = [ones(1,size(Pos,1)); Pos(:,1)'; Pos(:,2)'];
D = [zeros(3,3) D_];

c = zeros(1,3);

%Nonstruction von der Ausgleichung
N = [A'*A D';
     D  zeros(3,3)];
n = [A'*y; zeros(3,1)];

x_lam_dach = inv(N)*n;
x_dach = x_lam_dach(1:10)
lambda_dach = x_lam_dach(11:end);

%b
for i = 1:size(Pos,1)
    for k = 1:2
        delta_lambda= (Pos(i,1) - P_int(k,1))^2;
        delta_phi = (Pos(i,2) - P_int(k,2))^2;

        r = sqrt(delta_lambda + delta_phi);
        
        %festlegung für r = 0 
        if r == 0
            PHI_(i,k) = 0;
        else
            PHI_(i,k) = r^2 * log(r);
        end
    end
end

RBF = x_dach(4:end)'*PHI_;
h_int = x_dach(1) + x_dach(2)*P_int(:,1) + x_dach(3)*P_int(:,2) + RBF';

%% Aufgabe 2

 

       
 
 
 
 
 
 
 
 
 
 
 
 
 
 