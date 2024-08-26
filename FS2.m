% Figures S2 & S3.
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

% HF Map.
figure(53); clf; hold on;
plot(Blon,Blat,'-k');
plot(-68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
scatter(EQlon,EQlat,Rs/5,'o','MarkerEdgeColor','k','MarkerFaceColor',GREY);
scatter([S.Slon],[S.Slat],25,Thf,'filled','d','MarkerEdgeColor','b');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
xlabel('Longitude'); ylabel('Latitude'); 
colormap('turbo'); cb=colorbar(); clim([min(Thf) max(Thf)]);
ylabel(cb,'Time (year)');
pbaspect([1 Va 1]);

% EQ Map.
figure(52); clf; hold on;
plot(Blon,Blat,'-k');
plot(-68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
scatter([S.Slon],[S.Slat],10,'bd');
scatter(EQlon,EQlat,Rs,Teq,'o','filled','MarkerEdgeColor','k');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
xlabel('Longitude'); ylabel('Latitude'); 
colormap('turbo'); cb=colorbar(); clim([min(Thf) max(Thf)]);
ylabel(cb,'Time (year)');
pbaspect([1 Va 1]);





%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

