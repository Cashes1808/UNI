clc
clear all
close all
%% Einlesen del Bilder 
I_1=double(imread("cameraman.tif"));
I_2=double(imread("stuttgart_laser.tif"));
%% Aufgabe 1

% Boxoperator
I_1_n_B = falte(I_1,"Box");

%Kontrolle
A=(1/9)*ones(3,3);
I_new_box_kontrolle=uint8(imfilter(I_1,A));

% Visualisierung
figure
subplot(1,3,1)
I_1_show=imshow(uint8(I_1));
title("Originalbild")
subplot(1,3,2)
I_1_n_show=imshow(uint8(I_1_n_B));
title("Mit Box-Operator Gefaltetes Bild")
subplot(1,3,3)
imshow(I_new_box_kontrolle), title('Boxfilter Kontrolle')

% Sobel - x
I_1_n_B = falte(I_1,"Sobel-x");
%Kontrolle
B=(1/8)*[1 0 -1;2 0 -2;1 0 -1];
I_new_sobelx_kontrolle=uint8(imfilter(I_1,B));

% Visualisierung
figure
subplot(1,3,1)
I_1_show=imshow(uint8(I_1));
title("Originalbild")
subplot(1,3,2)
I_1_n_show=imshow(uint8(I_1_n_B));
title("Mit Sobel-x Gefaltetes Bild")
subplot(1,3,3)
imshow(I_new_sobelx_kontrolle), title('Sobel-x Kontrolle')

% Sobel - y
I_1_n_B = falte(I_1,"Sobel-y");

%Kontrolle
C=(1/8)*[1 2 1; 0 0 0; -1 -2 -1];
I_new_sobely_kontrolle=uint8(imfilter(I_1,C));

% Visualisierung
figure
subplot(1,3,1)
I_1_show=imshow(uint8(I_1));
title("Originalbild")
subplot(1,3,2)
I_1_n_show=imshow(uint8(I_1_n_B));
title("Mit Sobel-y Gefaltetes Bild")
subplot(1,3,3)
imshow(I_new_sobely_kontrolle), title('Sobel-y Kontrolle')




% Laplace
I_1_n_B = falte(I_1,"Laplace");
%Kontrolle
D=[0 1 0;1 -4 1;0 1 0];
I_new_laplace_kontrolle=uint8(imfilter(I_1,D));


% Visualisierung
figure
subplot(1,3,1)
I_1_show=imshow(uint8(I_1));
title("Originalbild")
subplot(1,3,2)
I_1_n_show=imshow(uint8(I_1_n_B));
title("Mit Laplace Gefaltetes Bild")
subplot(1,3,3)
imshow(I_new_laplace_kontrolle), title('Laplace Kontrolle')

%% Aufgabe 2
%Hoche neigungen bestimmen
I_2_sobelx=falte(I_2,"Sobel-x");
I_2_sobely=falte(I_2,"Sobel-y");

for i=1:size(I_2,1)
    for k=1:size(I_2,2)
        dZ_x(i,k)=I_2_sobelx(i,k);
        dZ_y(i,k)=I_2_sobely(i,k);
        N(i,k)=sqrt(dZ_y(i,k)^2+dZ_x(i,k)^2);
    end
end

N(N>4)=255;
N(N<4)=0;

figure
imshow(uint8(N))
title("Binärbild Höhenmodel")
%% funktionen
function [I_1_n] = falte(I_1,Filtername)
%Die Funktion "falte" faltet eine Matrix I_1 mit über der Funktion "falter".
% Dafür ist eine Liste an Auswahlmatrizen vordefiniert die über "Filtername" erreicht werden können.
%
%IN: 
%   I_1: n mal n Matrix das geflatet werden soll.
%   Filtername:
%               "Box": Wendet den Boxoperator (1/9)*[Ones(3x3)]
%               "Sobel-x": Wendet den Sobel Filter in x richung an. 
%               "Sobel-y": Wendet den Sobel Filter in y richung an. 
%   	        "Laplace": Wendet den Laplace Filter an. 
%
%
%OUT:
%   I_1_n:                 Die gefaltete Matrix.
%
%
%Author:                   Domonkos Veress
%Last update:              05.12.2022
%Project:                  Bildverarbeitung

% Anlage der vordefinierten Matrixen
%box-operator
if Filtername == "Box"
 f_streck_ecke=1/4; 
 Auswahlmatrix_ecke_1=[1 1;1 1];
 Auswahlmatrix_ecke_2=[1 1;1 1];
 Auswahlmatrix_ecke_3=[1 1;1 1];
 Auswahlmatrix_ecke_4=[1 1;1 1];

 f_streck_kante=1/6;
 Auswahlmatrix_Zeile_1=[1,1,1;1,1,1];
 Auswahlmatrix_Zeile_l=[1,1,1;1,1,1];
 Auswahlmatrix_Spalte_1=[1,1,1;1,1,1];
 Auswahlmatrix_spalte_l=[1,1,1;1,1,1];

f_streck_mitte=1/9; 
Auswahlmatrix_mitte=[1,1,1;1,1,1;1,1,1];

I_1_n = falter(I_1, f_streck_ecke, Auswahlmatrix_ecke_1, Auswahlmatrix_ecke_2, Auswahlmatrix_ecke_3, Auswahlmatrix_ecke_4, ...
                               f_streck_kante, Auswahlmatrix_Zeile_1, Auswahlmatrix_Zeile_l, Auswahlmatrix_Spalte_1, Auswahlmatrix_spalte_l, ...
                               f_streck_mitte, Auswahlmatrix_mitte);

end


%Sobel-operator in x-Richtung
if Filtername == "Sobel-x"
 f_streck_ecke=1/8; 
 Auswahlmatrix_ecke_4=[1 0;2 0];
 Auswahlmatrix_ecke_3=[0 -1;0 -2];
 Auswahlmatrix_ecke_2=[2 0;1 0];
 Auswahlmatrix_ecke_1=[0 -2;0 -1];

 f_streck_kante=1/8;
 Auswahlmatrix_Zeile_1=[1,0,-1;2,0,-2];
 Auswahlmatrix_Zeile_l=[2,0,-2;1,0,-1];
 Auswahlmatrix_Spalte_1=[1,0,2;0,1,0];
 Auswahlmatrix_spalte_l=[0,-1,0;-2,0,-1];

f_streck_mitte=1/8; 
Auswahlmatrix_mitte=[1,0,-1;2,0,-2;1,0,-1];

I_1_n = falter(I_1, f_streck_ecke, Auswahlmatrix_ecke_1, Auswahlmatrix_ecke_2, Auswahlmatrix_ecke_3, Auswahlmatrix_ecke_4, ...
                               f_streck_kante, Auswahlmatrix_Zeile_1, Auswahlmatrix_Zeile_l, Auswahlmatrix_Spalte_1, Auswahlmatrix_spalte_l, ...
                               f_streck_mitte, Auswahlmatrix_mitte);

end

%Sobel-operator in y-Richtung
if Filtername == "Sobel-y"
 f_streck_ecke=1/8; 
 Auswahlmatrix_ecke_4=[1 2;0 0];
 Auswahlmatrix_ecke_3=[2 1;0 0];
 Auswahlmatrix_ecke_2=[0 0;-1 -2];
 Auswahlmatrix_ecke_1=[0 0;-2 -1];

 f_streck_kante=1/8;
 Auswahlmatrix_Zeile_1=[1,2,1;0,0,0];
 Auswahlmatrix_Zeile_l=[0,0,0;-1,-2,-1];
 Auswahlmatrix_Spalte_1=[1,2,0;0,-1,-2];
 Auswahlmatrix_spalte_l=[2,1,0;0,-2,-1];

f_streck_mitte=1/8; 
Auswahlmatrix_mitte=[1,2,1;0,0,0;-1,-2,-1];

I_1_n = falter(I_1, f_streck_ecke, Auswahlmatrix_ecke_1, Auswahlmatrix_ecke_2, Auswahlmatrix_ecke_3, Auswahlmatrix_ecke_4, ...
                               f_streck_kante, Auswahlmatrix_Zeile_1, Auswahlmatrix_Zeile_l, Auswahlmatrix_Spalte_1, Auswahlmatrix_spalte_l, ...
                               f_streck_mitte, Auswahlmatrix_mitte);

end



%Laplace-operator
if Filtername == "Laplace"
 f_streck_ecke=1; 
 Auswahlmatrix_ecke_4=[0 1;1 -4];
 Auswahlmatrix_ecke_3=[1 0;-4 1];
 Auswahlmatrix_ecke_2=[1 -4;0 1];
 Auswahlmatrix_ecke_1=[-4 1;1 0];

 f_streck_kante=1;
 Auswahlmatrix_Zeile_1=[0,1,0;1,-4,1];
 Auswahlmatrix_Zeile_l=[1,-4,1;0,1,0];
 Auswahlmatrix_Spalte_1=[0,1,1;-4,0,1];
 Auswahlmatrix_spalte_l=[1,0,-4;1,1,0];

f_streck_mitte=1; 
Auswahlmatrix_mitte=[0,1,0;1,-4,1;0,1,0];

I_1_n = falter(I_1, f_streck_ecke, Auswahlmatrix_ecke_1, Auswahlmatrix_ecke_2, Auswahlmatrix_ecke_3, Auswahlmatrix_ecke_4, ...
                               f_streck_kante, Auswahlmatrix_Zeile_1, Auswahlmatrix_Zeile_l, Auswahlmatrix_Spalte_1, Auswahlmatrix_spalte_l, ...
                               f_streck_mitte, Auswahlmatrix_mitte);

end

end
function [I_1_n] = falter(I_1, f_streck_ecke, Auswahlmatrix_ecke_1, Auswahlmatrix_ecke_2, Auswahlmatrix_ecke_3, Auswahlmatrix_ecke_4, ...
                               f_streck_kante, Auswahlmatrix_Zeile_1, Auswahlmatrix_Zeile_l, Auswahlmatrix_Spalte_1, Auswahlmatrix_spalte_l, ...
                               f_streck_mitte, Auswahlmatrix_mitte)

%Die Funktion "falter" wendet auf eine Matrix I_1 eine auswahl an Matrizen
%uns Streckungsfaktoren an. 
%Die Faltung im inneren der Matrix ist unkopliziert. Eine Faltungsmaske kann ohne Probleme angewendet werden.
%Die faltung an den Rändern sowie den Ecken wird so gelöst, dass nur die umrandenden Elemente berücksichtigt werden.
%
%IN:
%   I_1: eine n mal n Matrix die gefaltete werden soll. 
%  
%   f_streck_ecke:        Der Streckungsfaktor an den Ecken. 
%   Auswahlmatrix_ecke_1: Auswahlmatrix an der linken oberen Ecke. Die
%                         Matrix ist [2x2] groß. 
%   Auswahlmatrix_ecke_2: Auswahlmatrix an der rechten oberen Ecke. Die
%                         Matrix ist [2x2] groß.
%   Auswahlmatrix_ecke_3: Auswahlmatrix an der linken unteren Ecke. Die
%                         Matrix ist [2x2] groß.
%   Auswahlmatrix_ecke_4: Auswahlmatrix an der rechten unteren Ecke. Die
%                         Matrix ist [2x2] groß.
%
%   
%   f_streck_kante:        Der Streckungsfaktor an den Kanten.
%   Auswahlmatrix_Zeile_1: Auswahlmatrix an der oberen Kante. Die
%                          Matrix ist [2x3] groß.
%   Auswahlmatrix_Zeile_l: Auswahlmatrix an der oberen Kante. Die
%                          Matrix ist [2x3] groß.
%   Auswahlmatrix_Spalte_1:Auswahlmatrix an der oberen Kante. Die
%                          Matrix ist [2x3] groß.
%   Auswahlmatrix_spalte_l:Auswahlmatrix an der oberen Kante. Die
%                          Matrix ist [2x3] groß.
%
%
%   f_streck_mitte:        Der Streckungsfaktor in der Mitte.
%   Auswahlmatrix_mitte:   Auswahlmatrix in der mitte der Matrix. Die
%                          Matrix ist [3x3] groß.
%
%
%OUT:
%   I_1_n:                 Die gefaltete Matrix.
%
%
%Author:                   Domonkos Veress
%Last update:              05.12.2022
%Project:                  Bildverarbeitung

%% Funktion 
%Vorbereitende definitionen 
d=size(I_1,1);
l=size(I_1,2);
I_1_n=zeros(d,l);

%berechnung der Ecken
I_1_n(1,1)=f_streck_ecke*(Auswahlmatrix_ecke_1(1,1)*I_1(1,1) + Auswahlmatrix_ecke_1(1,2)*I_1(1,2) + Auswahlmatrix_ecke_1(2,1)*I_1(2,1) + Auswahlmatrix_ecke_1(2,2)*I_1(2,2));
I_1_n(1,l)=f_streck_ecke*(Auswahlmatrix_ecke_2(1,1)*I_1(1,l-1) + Auswahlmatrix_ecke_2(1,2)*I_1(1,l) + Auswahlmatrix_ecke_2(2,1)*I_1(2,l-1) + Auswahlmatrix_ecke_2(2,2)*I_1(2,l));
I_1_n(d,1)=f_streck_ecke*(Auswahlmatrix_ecke_3(1,1)*I_1(d-1,1) + Auswahlmatrix_ecke_3(1,2)*I_1(d-1,2) + Auswahlmatrix_ecke_3(2,1)*I_1(d,1) + Auswahlmatrix_ecke_3(2,2)*I_1(d,2));
I_1_n(d,l)=f_streck_ecke*(Auswahlmatrix_ecke_4(1,1)*I_1(d-1,l-1) + Auswahlmatrix_ecke_4(1,2)*I_1(d-1,l) + Auswahlmatrix_ecke_4(2,1)*I_1(d,l-1) + Auswahlmatrix_ecke_4(2,2)*I_1(d,l));



%Berechnung der Kanten ohen den Ecken 
for k=2:l-1
    %Zeilen
    I_1_n(1,k)=f_streck_kante*(Auswahlmatrix_Zeile_1(1,1)*I_1(1,k-1) + Auswahlmatrix_Zeile_1(1,2)*I_1(1,k) + Auswahlmatrix_Zeile_1(1,3)*I_1(1,k+1) + Auswahlmatrix_Zeile_1(2,1)*I_1(2,k-1) + Auswahlmatrix_Zeile_1(2,2)*I_1(2,k) + Auswahlmatrix_Zeile_1(2,3)*I_1(2,k+1));
    I_1_n(d,k)=f_streck_kante*(Auswahlmatrix_Zeile_l(1,1)*I_1(d-1,k-1) + Auswahlmatrix_Zeile_l(1,2)*I_1(d-1,k) + Auswahlmatrix_Zeile_l(1,3)*I_1(d-1,k+1) + Auswahlmatrix_Zeile_l(2,1)*I_1(d,k-1) + Auswahlmatrix_Zeile_l(2,2)*I_1(d,k) + Auswahlmatrix_Zeile_l(2,3)*I_1(d,k+1));
end
for k=2:d-1
    %Zeilen
    I_1_n(k,1)=f_streck_kante*(Auswahlmatrix_Spalte_1(1,1)*I_1(k-1,1) + Auswahlmatrix_Spalte_1(1,2)*I_1(k-1,2) + Auswahlmatrix_Spalte_1(1,3)*I_1(k,1) + Auswahlmatrix_Spalte_1(2,1)*I_1(k,2) + Auswahlmatrix_Spalte_1(2,2)*I_1(k+1,1) + Auswahlmatrix_Spalte_1(2,3)*I_1(k+1,2));
    I_1_n(k,l)=f_streck_kante*(Auswahlmatrix_spalte_l(1,1)*I_1(k-1,l-1) + Auswahlmatrix_spalte_l(1,2)*I_1(k-1,l) + Auswahlmatrix_spalte_l(1,3)*I_1(k,l-1) + Auswahlmatrix_spalte_l(2,1)*I_1(k,l) + Auswahlmatrix_spalte_l(2,2)*I_1(k+1,l-1) + Auswahlmatrix_spalte_l(2,3)*I_1(k+1,l));
end


%berechnung der Mitte
for j=2:d-1
    for k=2:l-1
        I_1_n(j,k)=f_streck_mitte*(Auswahlmatrix_mitte(1,1)*I_1(j-1,k-1) + Auswahlmatrix_mitte(1,2)*I_1(j-1,k) + Auswahlmatrix_mitte(1,3)*I_1(j-1,k+1) + Auswahlmatrix_mitte(2,1)*I_1(j,k-1) + Auswahlmatrix_mitte(2,2)*I_1(j,k) + Auswahlmatrix_mitte(2,3)*I_1(j,k+1) + Auswahlmatrix_mitte(3,1)*I_1(j+1,k-1) + Auswahlmatrix_mitte(3,2)*I_1(j+1,k) + Auswahlmatrix_mitte(3,3)*I_1(j+1,k+1));
    end
end

end
