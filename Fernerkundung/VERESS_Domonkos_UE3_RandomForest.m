clc
clear all 
close all

%% Einstellungen 
Visuelle_Analysis = false;
onlyinclass = false;
%% Konstanten
    GlobalProgress_Zaehler = 12;
%% Progress bar einführen
%Aufgrund der Zeitintensität der aufgaben kann es sinn machen eine Progress
%bar einzuführen. 

GlobalProgress = waitbar(0/GlobalProgress_Zaehler,"Starting script");

%% Import data 
waitbar(1/GlobalProgress_Zaehler,GlobalProgress,"Importiere Daten")

%Bild laden
waitbar(2/GlobalProgress_Zaehler,GlobalProgress,"Bild einlesen: in progress")
I = double(imread("data\image.tif"));
%imshow(uint8(I(:,:,[3,1,11])))
waitbar(2/GlobalProgress_Zaehler,GlobalProgress,"Bild einlesen: Done")

%Info Laden
waitbar(3/GlobalProgress_Zaehler,GlobalProgress,"Info.mat laden: in progress")
I_info = load("data\info.mat");
waitbar(3/GlobalProgress_Zaehler,GlobalProgress,"Info.mat laden: done")

%Einlesen der Ground truth Daten
waitbar(4/GlobalProgress_Zaehler,GlobalProgress,"Ground truth laden: in progress")
GT = load("data\ground_truth.mat");
waitbar(4/GlobalProgress_Zaehler,GlobalProgress,"Ground truth laden: done")

%% Reorganisation der Daten
% In dieser sektion werden Daten so reorganisiert dass im weiter verlauf
% diese daten bearbeitet werden können. Die finale strucktur der
% reorganisierten daten sieht so aus : Pix_i = [kanal1:12, Klasse1:5]

%Vektorisierung des Bildes
waitbar(5/GlobalProgress_Zaehler,GlobalProgress,"Bild reorganisiert: in progress")
I_Neuorganisiert = Array_sizeMxNxK__2__Array_sizeMNxKx1(I);
waitbar(5/GlobalProgress_Zaehler,GlobalProgress,"Bild reorganisiert: Done")

%Vektorisierung der Ground Truth Daten
waitbar(6/GlobalProgress_Zaehler,GlobalProgress,"GT reorganisiert: in progress")
GT_Neuorganisiert = Array_sizeMxNxK__2__Array_sizeMNxKx1(GT.mask);
waitbar(6/GlobalProgress_Zaehler,GlobalProgress,"GT reorganisiert: Done")

%Zusammenfügen der einzelnen Datenquellen
waitbar(7/GlobalProgress_Zaehler,GlobalProgress,"DATA bilden: in progress")
DATA = [I_Neuorganisiert GT_Neuorganisiert];
waitbar(7/GlobalProgress_Zaehler,GlobalProgress,"DATA bilden: done")

CLASS1D = DATA(:,13:17)*[1:5]';

%% Untersuchungen vor dem Algorithmus

if Visuelle_Analysis == true
    Forest = GT.mask(:,:,1);
    %Jedes Layer der GT.mask matrix entpricht einer klasse
    % Nicht jedes Pixel kann zu eine der Klassen zugeordnet werden.
    


    %Eine untersuchung um das vorkommen einzelner gruawerte besser zu verstehen
    CH1 = DATA(:,1);
    CH1z1 = numel(find(CH1 == 1));
    CH1z57 = numel(find(CH1 == 57)); 
    CH1z112 = numel(find(CH1 == 112)); 
    CH1z187 = numel(find(CH1 == 187)); 
    CH1z210 = numel(find(CH1 == 210)); 
    % Schluss: Das anteilige Vorkommen der Faktoren (Grauwerte) im Kanal kann als Maß
    % für den Gewinn an Information dienen.
    % Der gefundene Wert als Schwellenwert.
end


if onlyinclass == true
    %% Hrasufiltern von Proben mit einem der 5 Kalssen
    waitbar(8/GlobalProgress_Zaehler,GlobalProgress,"Klassen filtern: in progress")
    class = isClass(DATA(:,13:17));
    cclass = find(class);
    waitbar(8/GlobalProgress_Zaehler,GlobalProgress,"Klassen filtern: done")
    %% Seprarating DATA in Klasse 
    waitbar(9/GlobalProgress_Zaehler,GlobalProgress,"DATA in Klasse: in progress")
    DATA_inclass = DATA(cclass,:);
    waitbar(9/GlobalProgress_Zaehler,GlobalProgress,"DATA in Klasse: done")
    
    %% DATA label in vektor
    waitbar(10/GlobalProgress_Zaehler,GlobalProgress,"DATA Label in Vektor: in progress")
    data_labeledvector = [];
    for i = 1:5
        rows = find(DATA_inclass(:,12 + i));
        class_name = ones(numel(rows),1)*i;
        DATA_inclass_labeled = [DATA_inclass(rows,1:12) class_name ];
        data_labeledvector = [data_labeledvector; DATA_inclass_labeled];
    end
    waitbar(10/GlobalProgress_Zaehler,GlobalProgress,"DATA Label in Vektor: done")
else 
    data_labeledvector = DATA;
end
%% TD und VD generieren
waitbar(11/GlobalProgress_Zaehler,GlobalProgress,"TD und VD generieren")

waitbar(11/GlobalProgress_Zaehler,GlobalProgress,"DATA2Random: in progress")
DATA_rndmzd = Matrix_randomizer(data_labeledvector);
waitbar(11/GlobalProgress_Zaehler,GlobalProgress,"DATA2Random: done")


waitbar(11/GlobalProgress_Zaehler,GlobalProgress,"Aufspaltung in TD und VD: in progress")
% decleration of training and validation data 
% I want to use 2/3 of GT for Training 

[row, ~, ~] = size(DATA_rndmzd);
n_TD = 0.5*floor(row/3);
%start_VD = n_TD + 1;

TD = DATA_rndmzd(1:n_TD,:);
VD = DATA_rndmzd(n_TD+1:row,:);
waitbar(11/GlobalProgress_Zaehler,GlobalProgress,"Aufspaltung in TD und VD: done")


%% Plot 

figure
plot(TD(:,1),TD(:,2), ".")
title ("Kanäle im featurespace")
hold on
plot(TD(:,1), TD(:,11), ".")
plot(TD(:,1), TD(:,4), ".")

legend('1/2', '1/11', '1/4')

hold off


if onlyinclass == true
    %% Generierung des Models 
    waitbar(12/GlobalProgress_Zaehler,GlobalProgress,"Model generieren: in progress")
     test = TreeBagger(100,TD(:,1:12),TD(:,13),...
                       OOBPredictorImportance="on", ...
                       Method="classification", ...
                       NumPredictorsToSample="all");
    waitbar(12/GlobalProgress_Zaehler,GlobalProgress,"Model generieren: done")
    
     %% Visualization of test
     
     view(test.Trees{1},Mode = "graph")
    
     figure
     plot(oobError(test))
     xlabel ("Gezogene Bäume")
     ylabel ("Out-of_Bag fehler")
    
    
     %% Validation
    waitbar(13/GlobalProgress_Zaehler,GlobalProgress,"Klassifizierung: in progress")
    val_test = predict(test,VD(:,1:12));
    val_array = str2num(cell2mat(val_test));
    
    Diff = VD(:,13)-val_array;
    n_of_errors = numel(find(Diff));
    n_of_errors_percent = (n_of_errors/numel(Diff))*100;
    
    
    
    %I_classfd = predict(test,I_Neuorganisiert);
    waitbar(13/GlobalProgress_Zaehler,GlobalProgress,"Klassifizierung: done")
end
%% Kalassifikation mit allen Kalssen

% Random Forest Classifier
Y_1D = double(TD(:,13:17)) * [1:5]';
Ntrees = 5;
BaggedEnsemble = TreeBagger(Ntrees, TD(:,1:12), Y_1D, 'OOBPred', 'On', 'Method', 'classification');

% Classification and results
VD_Class1D = VD(:,13:17)*[1:5]';
RFresult = predict(BaggedEnsemble, DATA(:,1:12));
RF = str2num(cell2mat(RFresult));

%% Berechnung konfusionsmatrix 
c_MAT = confusionmat(CLASS1D, RF);

for i = 1:6
    row_sum(i) = sum(c_MAT(i,:));
    col_sum(i) = sum(c_MAT(:,i));
end


%% Berechnung genauigkeit 
accuracy = sum(CLASS1D == RF) / numel(CLASS1D);

for i = 1:6
    Users_acc(i) = c_MAT(i,i)/row_sum(i);
    Produsers_acc(i) = c_MAT(i,i)/col_sum(i);
end

%% Visualisation

     view(BaggedEnsemble.Trees{1},Mode = "graph")
    
     figure
     plot(oobError(BaggedEnsemble))
     xlabel ("Gezogene Bäume")
     ylabel ("Out-of_Bag fehler")
     
     
     
     
 %% Bildbau 
[m, n, ~] = size(I);
I_class = reshape(RF,[m,n]);
     
figure
imshow(uint8(label2rgb(I_class)))

%% 
figure
imagesc(I_class)
%% ENDE
close(GlobalProgress)
%% Funktionen

function [array] = Matrixrectifier(matrix)
%This funktion takes the input matrix of size m x n and "rectifys it into
%an array of size m x 1. 
    [rows, ~] = size(matrix);
    
    array = [];
    for i = 1:rows
        array = [array; matrix(i,:)'];
    end
end

function [MATRIX_rectifyed_appended] = Array_sizeMxNxK__2__Array_sizeMNxKx1(Matrix, rangeofrectification)
%This funktion takes a Matrix of size n x m x k and creates a Matrix of
%size m x k x 1. 
%Optionally the number of 
    
    if ~exist("rangeofrectification","var")
        [~,~,depth] = size(Matrix);
        rangeofrectification = [1, depth];
    end
      

    MATRIX_rectifyed_appended = [];
    
        for i = rangeofrectification(1):rangeofrectification(2)
            temp_array = Matrixrectifier(Matrix(:,:,i));
            MATRIX_rectifyed_appended = [MATRIX_rectifyed_appended temp_array];
        end
end

function [final_class] = isClass(Data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
   
%progress bar 
    f = waitbar(0,"Classifying");
%Createsample
    total = Data + 1;
    n = 1;
    m = 10000;
    m_ = 10000;
    
    [row , col] = size(total);
    test = row/m;
    flo = floor(test);
    ttest = (test-flo)*m;
    final_class = [];
    for i = 1:flo
        
        sampler = total(n:m,:)*total(n:m,:)';
        sampel = diag(sampler);
    
        subtract_factor = col; %+ 3;
        
        preprefinal_class = sampel - subtract_factor;
    
        divisioning_factor = max(abs(preprefinal_class));

        if divisioning_factor == 0
            divisioning_factor = 1;
        end
        prefinal_class = preprefinal_class ./divisioning_factor;
        final_class = [final_class; abs(prefinal_class)];
        n = n+m_;
        m = m + m_;
    
        %Progress Bar Update
        waitbar(i/flo, f, sprintf('isClass done: %d %%', floor(i/flo*100)));
    end
    
    close (f)
end

function [rndm_Matrix] = Matrix_randomizer(Matrix)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    [row, ~] = size(Matrix);
    rndm = randperm(row);

    rndm_Matrix(1:row,:) = Matrix(rndm,:); 

end
