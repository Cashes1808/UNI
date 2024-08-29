clc
clear all
close all
%% Importieren der Daten
I1= double(imread("R0020851.jpg"));
I2= double(imread("R0020852.jpg"));

%Parameter I
I_C = 0.0216335; %[mm]
I_pixelgroesse = [0.0000055;0.0000055]; %[mm]
I_sensorgroesse = [4288; 2848]; %[pixel]
I_principal_point = [2143.5; 1423.5]; %[pixel]
I_intorr = [ -3933.36363636 0             2143.5 ;
             0              3933.36363636 1423.5 ;
             0              0             1     ];% bildkordinatensystem [K_matrix]

%Parameter I1

I1_extorr = [-0.827570749773   -0.557530411593      0.065471324026;
              0.549857636160	 -0.828569770064	 -0.105492730048;
       	      0.113062965097	 -0.051302790236	  0.992262460056]; % welt -> bildkoordinaten


I1_extorr_trans = [512980.995110	  5427701.526710    514.794290]';% Weltkoordinaten


%Parameter I2
I2_extorr = [-0.866544118547	 -0.490057510421	 -0.094577624687;
	          0.490759244246	 -0.871123481449	  0.017298677877;
	         -0.090866136699	 -0.031424776041	  0.995367182829]; % welt -> bildkoordinaten

I2_extorr_trans = [ 512966.499230	  5427725.063170	      517.490680]'; % Weltkoordinaten



%% Bildpunktmessung

% Messen der identschen Punkte in beiden Bildern
figure;
imshow(uint8(I1));
title('Punkt messen! (Bild 1)')
[x1,y1] = getpts;  % Misst Punkt [pixel]
close;

imshow(uint8(I2));
title('Punkt messen! (Bild 2)')
[x2,y2] = getpts;  % Misst Punkt [pixel]
close;



%% Räumlicher vorwärtschnitt Vektoriell

%Umwandlung in homogene Koordinaten
x1_pix = [x1 ; y1 ; 1];
x2_pix = [x2 ; y2 ; 1];

%Pixelkoordinaten in Kamerakoordinaten
x1_cam = inv(I_intorr) * x1_pix;
x2_cam = inv(I_intorr) * x2_pix;

V1 = I1_extorr' * x1_cam;
V2 = I2_extorr' * x2_cam; 

%berechnen von d, b, lambda
d = cross(V1,V2);
b = I2_extorr_trans - I1_extorr_trans;
lambda = [V1 V2 d]\b;

% berechnung der punktkoordinaten
P_0= I1_extorr_trans + lambda(1) * V1 + (1/2) * lambda(3) * d;
guete= lambda(3) * norm(d)

%% Räumlicher vorwärtschritt Projektiv

%Projektionsmatrix
P1= I_intorr * [I1_extorr,  -I1_extorr * I1_extorr_trans];
P2= I_intorr * [I2_extorr,  -I2_extorr * I2_extorr_trans];


A=[x1 * P1(3,:) - P1(1,:);
   y1 * P1(3,:) - P1(2,:);
   x2 * P2(3,:) - P2(1,:);
   y2 * P2(3,:) - P2(2,:)];

[U,D,V]= svd(A,0);
XX = V(:,end);
XYZ = XX / XX(4);
XYZ = XYZ(1:3);

%Vergeich der zwei methoden 

Diff = XYZ - P_0;


%% Projektion der signalized_points

%signalpunkte
p15 = [512868.940 5427723.833 280.885]';
p24 = [512982.899 5427801.281 323.376]';
p25 = [513047.601 5427704.443 326.998]';
p32 = [512781.152 5427721.195 229.286]';
p37 = [512866.691 5427556.801 229.433]';

%Gemessen projiziert
x_mess_bild1 = P1 * [XYZ;1];
x_mess_bild1 = x_mess_bild1/x_mess_bild1(3);

x_mess_bild2 = P2 * [XYZ;1];
x_mess_bild2 = x_mess_bild2/x_mess_bild2(3);

Diff_proj_mess = [x1,y1,x2,y2]' - [x_mess_bild1(1), x_mess_bild1(2), x_mess_bild2(1), x_mess_bild2(2)]'




% Transformation 
x_bild1 = P1 * [p24;1];
x_bild1 = x_bild1/x_bild1(3);
x_bild1 = x_bild1(1:2);

x_bild2= P2 * [p24;1];
x_bild2 = x_bild2/x_bild2(3);
x_bild2 = x_bild2(1:2);


%plot
figure
imshow(uint8(I1))
hold on
scatter([x_mess_bild1(1) , x_bild1(1)],[x_mess_bild1(2) , x_bild1(2)],150,'filled' )
title('Bild 1')
hold off

figure
imshow(uint8(I2))
hold on
scatter([x_mess_bild2(1) , x_bild2(1)],[x_mess_bild2(2) , x_bild2(2)],150,'filled' )
title('bild 2')
hold off



%%
A=[1 1  1;
   1  2 4;
   1 3 9];

y= [2,1,1]'

X_dach= inv(A'*A)*A'*y

