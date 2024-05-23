% Figure S1.
clear;

% Predefine some stuff.
welFILE1='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';
welFILE2='/Users/rschultz/Desktop/VacaMuerta/data/wells2.txt';

% Load the HF data.
load(welFILE1,'S');
W2=load(welFILE2);

% Load the Neuquén Basin boundaries.
data=load('../data/shapes_vaca_muerta/Neuquen.txt');
Blat=data(:,2); Blon=data(:,1);

% Plot.
%latB=[-40.5 -34.2]; lonB=[-72.0 -66.5]; % Basin.
latB=[-39.2 -37.3]; lonB=[-69.75 -67.55]; % HF wells.
Va=(max(latB)-min(latB))/(max(lonB)-min(lonB))/cosd(mean(latB));

% Map.
figure(51); clf; hold on;
plot(Blon,Blat,'-k');
plot(-68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
plot(W2(:,2),W2(:,3),'bd');
plot([S.Slon],[S.Slat],'bd','MarkerFaceColor','c');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
pbaspect([1 Va 1]);
