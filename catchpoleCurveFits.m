% Housekeeping
clear
clc

% Load files
L05 = load('Catchpole/0.5.csv');
C05 = polyfit(L05(:,1),L05(:,2),15);

L06 = load('Catchpole/0.6.csv');
C06 = polyfit(L06(:,1),L06(:,2),15);

L07 = load('Catchpole/0.7.csv');
C07 = polyfit(L07(:,1),L07(:,2),15);

L08 = load('Catchpole/0.8.csv');
C08 = polyfit(L08(:,1),L08(:,2),15);

L09 = load('Catchpole/0.9.csv');
C09 = polyfit(L09(:,1),L09(:,2),15);

L10 = load('Catchpole/1.0.csv');
C10 = polyfit(L10(:,1),L10(:,2),15);

L125 = load('Catchpole/1.25.csv');
C125 = polyfit(L125(:,1),L125(:,2),15);

L150 = load('Catchpole/1.5.csv');
C150 = polyfit(L150(:,1),L150(:,2),15);

L20 = load('Catchpole/2.0.csv');
C20 = polyfit(L20(:,1),L20(:,2),15);

%% Plots to chek validity
% disc = L05(:,1);
% 
% % 0.5
% plot(L05(:,1),L05(:,2))
% hold on
% plot (disc,polyval(C05,disc))
% 
% % 0.6
% plot(L06(:,1),L06(:,2))
% plot (disc,polyval(C06,disc))
% 
% % 0.7
% plot(L07(:,1),L07(:,2))
% plot (disc,polyval(C07,disc))
% 
% % 0.8
% plot(L08(:,1),L08(:,2))
% plot (disc,polyval(C08,disc))
% 
% % 0.9
% plot(L09(:,1),L09(:,2))
% plot (disc,polyval(C09,disc))
% 
% % 1.0
% plot(L10(:,1),L10(:,2))
% plot (disc,polyval(C10,disc))
% 
% % 1.25
% plot(L125(:,1),L125(:,2))
% plot (disc,polyval(C125,disc))
% 
% % 1.5
% plot(L150(:,1),L150(:,2))
% plot (disc,polyval(C150,disc))
% 
% % 2.0
% plot(L20(:,1),L20(:,2))
% plot (disc,polyval(C20,disc))

%% Save coefficients to .mat
save('catchpoleCoeffs.mat','C05','C06','C07','C08','C09','C10','C125','C150','C20');