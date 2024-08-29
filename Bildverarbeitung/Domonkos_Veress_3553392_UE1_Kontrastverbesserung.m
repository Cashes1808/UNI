clc
clear all
close all 

%% Einlesen des Bildes 
I=imread('pout.tif');
figure(1);imshow(I);title('Originalbild'); 

%% Ursprüngliches Histogramm
%Visualisierung
figure(2);imhist(I);title('Histogramm Originalbild'); 

%% Linear skaliertes Histogramm 
% Bestimmung max. und min. Gruawerte skallierung und offset
min_I=double (min(I(:)));
max_I=double (max(I(:)));
scale=( 255.0 - 0.)/( max_I - min_I);
offset=(max_I * 0.-min_I*255.)/(max_I - min_I);

% int8 zwischen 0 und 255 -> umformen in int16 (-32768 to 32767)
I_stretch=int16(I)* scale + offset;
% Rückwandlung von int16 in int8 
I_stretch = uint8 (I_stretch);

%Darstellung des Skalierten Bildes  
figure(3);imshow(I_stretch);title('Lineare Skalierung');
figure(4);imhist(I_stretch);title('Histogramm Lineare Skalierung');
% 
%% Berechnung Normiertes kummuliertes Histogramm

%kummuliertes, normiertes Histogramm
cumHist = cumsum(imhist(I)/numel(I));
% Visualisierung kummulatives Histogramm
x = linspace(0, 255, 256);
figure(5), plot(x,cumHist);title('kummulatives, normiertes Histogramm');
%Thresholds für die dunkelsten und hellsten Werte
thresh_low = 0.05; thresh_up = 1 - thresh_low;
% Bestimmung Position im Vektor cumHistogramm
min_I_t_h=max(find(cumHist <= thresh_low));
max_I_t_h=min(find(cumHist >= thresh_up));
 
 
%Widerhohlung der Skalierung 
% Bestimmung skallierung und offset Gruawerte
scale=( 255.0 - 0.)/( max_I_t_h - min_I_t_h);
offset=(max_I_t_h * 0.-min_I_t_h*255.)/(max_I_t_h - min_I_t_h);

% int8 zwischen 0 und 255 -> umformen in int16 (-32768 to 32767)
I_stretch=int16(I)* scale + offset;
% Rückwandlung von int16 in int8 
I_stretch = uint8 (I_stretch);

%Darstellung des Skalierten Bildes  
figure(6);imshow(I_stretch);title('Normiertes kummuliertes Histogramm');
figure(7);imhist(I_stretch);title('Histogramm Normiertes kummuliertes Histogramm');


%% Histogrammverebnung
I_neu_man=zeros(size(I),'double')
for i=1:256
    %Berechnung der verteilung
    count(i)=length(find(I==i));
    norm(i)=sum(count(1:i))/numel(I);
    k=find(I==i);

    %Neues Bild
    I_neu_man(k)=round(norm(i)*(numel(I)/256));
    I_neu_man_uint8=uint8(I_neu_man);
end

I_neu_man=I_neu_man(:)';

%Bilddarstellung nach der Verebnung
figure;imshow(I_neu_man_uint8);title('Verebnetes Bild');

%Erstellung des Histogramms
bins=256;
minI=min(I_neu_man); maxI=max(I_neu_man);
bin_width=(maxI-minI)/bins;
bin_num=1+floor((I_neu_man-minI)./bin_width);

for i=1:bins
    bincount(i)=nnz(bin_num==i);
end

%Vizualisierung des Histogramm
x=min_I+(0.5:bins-0.5).*bin_width;
figure;
bar(x,bincount)
title('histogrammverebnung')




