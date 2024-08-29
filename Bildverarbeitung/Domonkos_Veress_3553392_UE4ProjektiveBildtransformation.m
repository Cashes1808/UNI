pl,okz vclc
clear all 
close all
%% Einlesen der Bilder
I1_1=double(imread("Gebaeude_0004_half.jpg"));
I1_2=double(imread("Gebaeude_0005_half.jpg"));
I2_1=double(imread("R0020849.jpg"));
I2_2=double(imread("R0020850.jpg"));
I3_1=double(imread("Seminarraum1.jpg"));
I3_2=double(imread("Seminarraum2.jpg"));

%% Bildpaar 1

%Einlesen der Berechnungspunkte
[x_s,y_s, x_d,y_d,l_1] = Punktbestimmer(I1_1,I1_2);

%berechnung der Homographien
H_QuellZiel=H_Rechner(x_s,y_s,x_d,y_d);
H_ZielQuell=H_Rechner(x_d,y_d,x_s,y_s);

%interpolation in richtung Ziel
IN=H_interpolator(I1_1,H_QuellZiel);
figure
imshow(uint8(IN))


Interpolation in richtung Quelle
INN=H_interpolator(IN,H_ZielQuell);
figure
imshow(uint8(INN))


[INNN] = Panorama_kreator(I1_2,IN);
imshow(uint8(INNN))
%% Bildpaar 2

%Einlesen der Berechnungspunkte
[x_s,y_s, x_d,y_d,l_2] = Punktbestimmer(I2_1,I2_2);

%berechnung der Homographien
H_QuellZiel=H_Rechner(x_s,y_s,x_d,y_d);
H_ZielQuell=H_Rechner(x_d,y_d,x_s,y_s);

%interpolation in richtung Ziel
IN=H_interpolator(I2_1,H_QuellZiel);
figure
imshow(uint8(IN))


%Interpolation in richtung Quelle
INN=H_interpolator(IN,H_ZielQuell);
figure
imshow(uint8(INN))


[INNN] = Panorama_kreator(I2_2,INN);
imshow(uint8(INNN))

%% Bildpaar 3

%Einlesen der Berechnungspunkte
[x_s,y_s, x_d,y_d,l_3] = Punktbestimmer(I3_1,I3_2);

%berechnung der Homographien
H_QuellZiel=H_Rechner(x_s,y_s,x_d,y_d);
H_ZielQuell=H_Rechner(x_d,y_d,x_s,y_s);

%interpolation in richtung Ziel
IN=H_interpolator(I3_1,H_QuellZiel);
figure
imshow(uint8(IN))


%Interpolation in richtung Quelle
INN=H_interpolator(IN,H_ZielQuell);
figure
imshow(uint8(INN))


[INNN] = Panorama_kreator(I3_2,INN);
imshow(uint8(INNN))

%% Funktionen

function [xi,yi] = Einleser(I_toread)
% Die funktion einleser ermöglicht es koordinaten von Pixel von einem bild
% mit der Maus einzulesen. 
%   IN: 
%       I_toread: Das bild auf welchem Koordinaten betsimmt werden sollen.
%
%   OUT: 
%       xi: Vektor der x-Koordinaten der eingelesenen Punkte
%       yi: Vektor der Y-Koordinaten der eingelesenen Punkte
%   Author: Domonkos Veress 
%   Last Update: 05.01.2023


I_source_figure=imtool(uint8(I_toread));
[xi,yi]=getline(I_source_figure);
close(I_source_figure)

end

function [H] = H_Rechner(x_s, y_s, x_d, y_d)

%Die Funktion H_Rechner berechent eine Homographie von einem Quellbild
%in ein Zielbild.
%IN:
%   x_s: Vekor der x-Koordinaten des Quellbildes
%   y_s: Vekor der y-Koordinaten des Quellbildes
%   x_d: Vekor der x-Koordinaten des Zielbildes
%   y_d: Vekor der y-Koordinaten des Zielbildes
%OUT:
%   H: Matrix der Homographie
%Author: Domonkos Veress
%Last updated: 05.01.2023


% k=1;
% for i=1:length(x_s)
% A(k:k+1,:)=[x_s(i)  y_s(i)  1  0       0       0  -x_d(i)*x_s(i)  -x_d(i)*y_s(i)  -x_d(i);
%             0        0      0  x_s(i)  y_s(i)  1  -y_d(i)*x_s(i)  -y_d(i)*y_s(i)  -y_d(i)];
% k=k+2;
% end

k=1;
for i=1:length(x_s)
A(k:k+1,:)=[0        0      0  x_s(i)  y_s(i)  1  -y_d(i)*x_s(i)  -y_d(i)*y_s(i)  -y_d(i)
            x_s(i)  y_s(i)  1  0       0       0  -x_d(i)*x_s(i)  -x_d(i)*y_s(i)  -x_d(i)];
k=k+2;
end
[Eigenvektoren,Eigenwerte]=eig(A'*A);
[~,ind]=sort(diag(Eigenwerte));
Eigenvektoren=Eigenvektoren(:,ind);
H_h=Eigenvektoren(:,1);
H=[H_h(1) H_h(2) H_h(3);
   H_h(4) H_h(5) H_h(6);
   H_h(7) H_h(8) H_h(9)];
end

function [H_x_d] = H_Transformator(x_s, y_s,H)

%Die funktion H_Transformator berechnet eine Transformation anhand der
%Matrix H und den Vektoren x_s und y_s, diese werden von der Funktion zu Homogenen Koordinaten umgewandelt.
%IN:
%   x_s: Vekor der x-Koordinaten des Quellbildes
%   y_s: Vekor der y-Koordinaten des Quellbildes
%   H: Transformationsmatrix
%OUT: 
%   H_x_d: Die transformierten Koordinaten in einer Matrix. Die Spalten
%   entsprechen den einzelenen Koordinatenpaaren. 
%Author: Domonkos Veress
%Last updated: 05.01.2023

for i=1:length(x_s)
    p_s=[x_s(i); y_s(i); 1];
    p_d_transformiert=H*p_s;
    H_x_d(:,i)=p_d_transformiert/p_d_transformiert(3);
end

end

function [x_s,y_s,x_d,y_d,l] = Punktbestimmer(I_Quelle,I_Ziel)

%Die Funktion Punktbestimmer ermöglicht das einlesen von Pixeln von einem
%Quellbild und einem Zielbild sowie das Anzeigen der Eingelesenen Punkte
%und der homographisch Transformierten Punkte im Zielbild. 
%IN: 
%   I_Quelle: Quellbild
%   I_Ziel: Zielbild
%OUT: 
%   x_s: Vekor der x-Koordinaten des Quellbildes
%   y_s: Vekor der y-Koordinaten des Quellbildes
%   x_d: Vekor der x-Koordinaten des Zielbildes
%   y_d: Vekor der y-Koordinaten des Zielbildes
%   l: Mass der sicherheit.
%Author: Domonkos Veress
%Last updated: 05.01.2023



%I_Quelle
[x_s,y_s] = Einleser(I_Quelle);

%I_Ziel
[x_d,y_d] = Einleser(I_Ziel);

%funktion um H zu berechnen
H=H_Rechner(x_s,y_s,x_d,y_d);
%funktion um H*x_d zu berechnen
H_x_d=H_Transformator(x_s,y_s,H);
l=[x_d-H_x_d(1,:)';y_d-H_x_d(2,:)'];

%Vizualization
figure;
imshow(uint8(I_Quelle))
hold on
plot(x_s,y_s,'r+','MarkerSize',30)
hold off

figure;
imshow(uint8(I_Ziel))
hold on
plot(x_d,y_d,'r+','MarkerSize',30)
plot(H_x_d(1,:),H_x_d(2,:),'b+','MarkerSize',30)
legend('Eingegebene Punkte','Transformierte Punkte von Quell nach Ziel')
hold off
end

function [I_interpoliert] =H_interpolator(I_Quelle,H)

%Die funktion H_interpolator interpolier farbwerte nach der H_transformation 
%IN: 
%   I_quelle: Das Bild welches transformiert werden soll.
%   H: Transformationsmatrix
%OUT:
%   I_interpoliert: Das bild welches nach der transformation interpoliert wurde.
%Author: Domonkos Veress
%Last updated: 05.01.2023

[width, height, ~] = size(I_Quelle);

[xi,yi] = meshgrid(1:height, 1:width);

TransPoints = [xi(:) yi(:) ones(length(yi(:)),1)]';

I_transformiert_normiert=H_Transformator(TransPoints(1,:), TransPoints(2,:),H);

xi = reshape(I_transformiert_normiert(1,:),width,height);
yi = reshape(I_transformiert_normiert(2,:),width,height);
I_interpoliert(:,:,1) = interp2(I_Quelle(:,:,1),xi,yi,'linear',0); 
I_interpoliert(:,:,2) = interp2(I_Quelle(:,:,2),xi,yi,'linear',0);
I_interpoliert(:,:,3) = interp2(I_Quelle(:,:,3),xi,yi,'linear',0);
end

function [I_Panorama] = Panorama_kreator(I_in1,I_in2)
%Die funktion Panorama_kreator vereinigt Bilder im um ein Panoramabild zu
%berechnen.
%IN:
%   I_in1: das bild auf welches die Panoramarechnung gemacht werden soll.
%   I_in2: Das bild welches auf das Referenzbild gerechnet werden soll
%Author: Domonkos Veress
%Last updated: 05.01.2023

%Addition zum Panoramabild
I_Panorama(:,:,1)=I_in1(:,:,1)+I_in2(:,:,1)/2;
I_Panorama(:,:,2)=I_in1(:,:,2)+I_in2(:,:,2)/2;
I_Panorama(:,:,3)=I_in1(:,:,3)+I_in2(:,:,3)/2;
end