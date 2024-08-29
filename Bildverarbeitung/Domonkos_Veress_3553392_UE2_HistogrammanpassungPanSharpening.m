clc 
clear all
close all

format long g
%% Laden der Bilder
I_1=imread("Histogrammanpassung_img_1.jpg");
I_2=imread("Histogrammanpassung_img_2.jpg");

%% Histogramanpassung
% Das Histogram von I_2 soll auf das Histogram von I_1 angepasst werden. 
%Dazu muss erstmal das komulative Histogram von Bild 1 sowie Bild 2 erzeugt
%werden.

%berechnen des Bildhistograms und des Kommulierten Hitograms von Bild 1
%[f(I^1)]
I_1_hist=imhist(I_1);
I_1_normcumhist=cumsum(I_1_hist/numel(I_1));


%berechnen des Bildhistograms und des Kommulierten Hitograms von Bild 2
%[g(I^1)]
I_2_hist=imhist(I_2);
I_2_normcumhist=cumsum(I_2_hist/numel(I_2));

%Test
numel_I_1=numel(I_1);
numel_I_1_hist=sum(I_1_hist);
numel_I_1_cumhist=I_1_normcumhist(end);

numel_I_2=numel(I_2);
numel_I_2_hist=sum(I_2_hist);
numel_I_2_cumhist=I_2_normcumhist(end);

if numel_I_1_hist || numel_I_1_cumhist ~= numel(I_1)
elseif numel_I_2_hist ||numel_I_2_cumhist ~= numel(I_2)
    disp("iwas ist falsch im Hiastogram")
end


%Berechnung der ...
 %Suchen aller Punkte in Bild 2 mit bestimmtem Grauwert i
 %Bestimmung der kumulierten Wahrscheinlichkeit aller Grauwerte <=i in Bild 2
 %Bestimmung der Grauwerte mit ähnlichster kumulierter Wahrscheinlichkeit in Bild 1
 %Zuordnung der neuen Grauwerte in Bild 2
I_2_n_help=zeros(size(I_2));
 for i=0:255
    k=find(I_2==i);                      
    v=I_2_normcumhist(i+1);                
    [~,e]=min(abs(v-I_1_normcumhist));     
    I_2_n_help(k)=e-1;                       
end

%Neues I_2
I_2_n=uint8(I_2_n_help);

%Test der Histogramanpassung
I_2_histeq=histeq(I_2,imhist(I_1));

%Visualisierung
figure;
subplot(1,2,1)
imshow(I_1)
title ("Bild I_1 ohne Anpassung")
subplot(1,2,2)
imhist(I_1)
title ("Histogramm I_1 ohne Anpassung")


figure;
subplot(1,2,1)
imshow(I_2)
title ("Bild I_2 ohne Anpassung")
subplot(1,2,2)
imhist(I_2)
title ("Histogramm I_2 ohne Anpassung")


figure;
subplot(1,2,1)
imshow(I_2_n)
title("I_2 nach Anpassung")
subplot(1,2,2)
imhist(I_2_n)
title(" I_2 Histogramm nach Anpassung")
 
 
figure;
subplot(1,2,1)
imshow(I_2_histeq)
title("I_2 nach histeq")
subplot(1,2,2)
imhist(I_2_histeq)
title("Histogramm nach histeq")
%% Histogramverebnung im farbraum
I_3=double(imread('Histogrammanpassung_img_1_RGB.jpg'));

%Normierung
R=I_3(:,:,1)./255;
G=I_3(:,:,2)./255;
B=I_3(:,:,3)./255;


%hue
%H=acosd(1./2*(R-G+R-B)./sqrt(((R-G).^2+(R-B).*(G-B))));
lol=1/2.*(R-G+R-B)./sqrt(((R-G).^2+(R-B).*(G-B)));
lol(lol>1)=1;  
H=acosd(lol); 
H(B>G) = 360-H(B>G);
H(isnan(H))=0;

%saturation
S=1-((3./(R+G+B)).*min(B,min(R,G)));   
S(isnan(S))=0;

%intensity
I=(R+B+G)./3;

%Histogramverebnung
I_verebnet=zeros(size(I));
I_100=uint8(I*255);
norm_cum_hist_I=cumsum(imhist(I_100)/numel(I_100));

for g=1:255
k=find(I_100 == g);
I_verebnet(k) = norm_cum_hist_I(g);
end

I_verebnet100=I_verebnet*255;
I_verebnet100_u8=uint8(I_verebnet100);


%Rücktransformation in RGB-Raum 
I_einsatz=I_verebnet;


R_n=zeros(size(R));
G_n=zeros(size(G));
B_n=zeros(size(B));

% für H<=120°:

i=(H<=120);	
B_n(i)=I_einsatz(i).*(1-S(i));
R_n(i)=I_einsatz(i).*(1+S(i).*cosd(H(i))./cosd(60-H(i)));
G_n(i)=3.*I_einsatz(i)-(R_n(i)+B_n(i));


%für 120°<H<=240°:

i=(120<H & H<=240);
H(i)=H(i)-120;
R_n(i)=I_einsatz(i).*(1-S(i));
G_n(i)=I_einsatz(i).*(1+S(i).*cosd(H(i))./cosd(60-H(i)));
B_n(i)=3.*I_einsatz(i)-(R_n(i)+G_n(i));


%für H>240°:

i=(240<H);
H(i)=H(i)-240;
G_n(i)=I_einsatz(i).*(1-S(i));
B_n(i)=I_einsatz(i).*(1+(S(i).*cosd(H(i)))./cosd(60-H(i)));
R_n(i)=3.*I_einsatz(i)-(G_n(i)+B_n(i));

%Zusammensetzen des neuen Bildes:
I_n=zeros(size(I_3));
I_n(:,:,1)=255*R_n;
I_n(:,:,2)=255*G_n;
I_n(:,:,3)=255*B_n;
I_n_u8=uint8(I_n);


figure;
subplot(2,2,1)
imshow(uint8(I_3))
title("original RGB")
subplot(2,2,2)
imhist(I_100), 
title('Intensity Orginal')
subplot(2,2,3)
imshow(I_n_u8); 
title("Verebnetes RGB Bild")
subplot(2,2,4)
imhist(I_verebnet100_u8),
title('Intensity verebnet')



%% Pan-Sharpening

I_4_lowRGB = imread('DMC_haeuser_lowRGB_original.bmp');
I_4_bw = double(imread('DMC_haeuser_pan.bmp'));

%Anpassen der Farbbildgröße an das des panchromatischen Bildes
I_4_lowRGB=im2double(imresize(I_4_lowRGB,[396,395]));

R=I_4_lowRGB(:,:,1);		
G=I_4_lowRGB(:,:,2);		
B=I_4_lowRGB(:,:,3);		

%Konvertierung in HSI
H=acosd(1./2*(R-G+R-B)./sqrt(((R-G).^2+(R-B).*(G-B)))); 
H(B>G) = 360-H(B>G);	
H(isnan(H))=0; 

S=1-3./(R+G+B).*min(min(R,G),B);
S(isnan(S))=0;	

I=(R+G+B)/3;

R_n=zeros(size(R));
G_n=zeros(size(G));
B_n=zeros(size(B));

%Rücktransformation in RGB-Raum für H<=120°:

i=(H<=120);	
B_n(i)=I_4_bw(i).*(1-S(i));
R_n(i)=I_4_bw(i).*(1+S(i).*cosd(H(i))./cosd(60-H(i)));
G_n(i)=3.*I_4_bw(i)-(R_n(i)+B_n(i));


%Rücktransformation in RGB-Raum für 120°<H<=240°:

i=(120<H & H<=240);
H(i)=H(i)-120;
R_n(i)=I_4_bw(i).*(1-S(i));
G_n(i)=I_4_bw(i).*(1+S(i).*cosd(H(i))./cosd(60-H(i)));
B_n(i)=3.*I_4_bw(i)-(R_n(i)+G_n(i));


%Rücktransformation in RGB-Raum für H>240°:
i=(240<H);
H(i)=H(i)-240;
G_n(i)=I_4_bw(i).*(1-S(i));
B_n(i)=I_4_bw(i).*(1+(S(i).*cosd(H(i)))./cosd(60-H(i)));
R_n(i)=3.*I_4_bw(i)-(G_n(i)+B_n(i));

%Zusammensetzen des neuen Bildes:
I_n=zeros(size(I_4_bw));
I_n(:,:,1)=R_n;
I_n(:,:,2)=G_n;
I_n(:,:,3)=B_n;
I_n=uint8(I_n);

figure;
subplot(2,2,1)
imshow(uint8(I_4_bw))
title ("Pnachromatisches Bild I_4 ohne Anpassung")
subplot(2,2,3)
imshow(I_4_lowRGB)
title ("RGB Bild I_4 ohne Anpassung")
subplot(2,2,[2,4])
imshow(I_n)
title("Bild nach pan sharpening")

