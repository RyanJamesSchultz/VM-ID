% Figures S1.
clear;

% Predefine some stuff.
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/CAT.mat';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';

% Load the EQ catalogue & HF data.
load(catFILE,'EQlat','EQlon','EQdep','EQtime','EQmag');
load(welFILE,'S');

% Load the Neuquén Basin boundaries.
data=load('../data/shapes_vaca_muerta/Neuquen.txt');
Blat=data(:,2); Blon=data(:,1);

% Get decimal years.
Teq=years(EQtime-datetime(2016,01,01))+2016;
Thf=years(arrayfun(@(s) mean(s.T),S)-datetime(2016,01,01))+2016;

% Make the magnitude size axis.
Rs=getMscale(EQmag)/10;

% Plotting setup.
%latB=[-40.5 -34.2]; lonB=[-72.0 -66.5]; % Basin.
latB=[-39.2 -37.3]; lonB=[-69.75 -67.55]; % HF wells.
Va=(max(latB)-min(latB))/(max(lonB)-min(lonB))/cosd(mean(latB));
GREY=[0.85,0.85,0.85];

% Station coordinates.
latS1=[-38.21 -38.38 -38.60 -38.28 -38.30];
lonS1=[-69.01 -69.13 -69.10 -68.47 -68.65];
latS2=[-38.40 -38.54 -38.44 -38.81 -38.69 -38.30];
lonS2=[-68.64 -68.54 -68.92 -68.85 -68.50 -68.76];
latS3=[-38.3477  -37.8875 -36.5877];
lonS3=[-68.7871  -71.0620 -69.0653];

% EQ Map.
figure(51); clf; hold on;
plot(Blon,Blat,'-k');
plot(-68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
scatter([S.Slon],[S.Slat],5,'bd');
scatter(EQlon,EQlat,Rs,'o','MarkerEdgeColor','k','MarkerFaceColor',GREY);
plot(lonS1,latS1,'^','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','#1DF502');
plot(lonS2,latS2,'v','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','#B1FAA7');
plot(lonS3,latS3,'>','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','#2FA120');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
xlabel('Longitude'); ylabel('Latitude'); 
%colormap('turbo'); cb=colorbar(); clim([min(Thf) max(Thf)]);
%ylabel(cb,'Time (year)');
pbaspect([1 Va 1]);





%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

