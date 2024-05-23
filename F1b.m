% Script to draw in parts of Figure 1b.
clear;

% Predefine some stuff.
dM=0.3;
Mc=2.5;
dR=[10.0 20.0];
dTs=30;
Scut=1;
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/20240314_VM_quakes_all_soruces.csv';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';

% Boundary polygons and stuff.
latB=[-39.2 -37.3 -37.3 -39.2];
lonB=[-69.8 -69.8 -66.8 -66.8];
depB=-1;
tB=-1;
mB=-1;

% Load the EQ catalogue data.
[EQlat,EQlon,EQdep,EQtime,EQmag]=parseINPRES(catFILE,latB,lonB,depB,tB,mB);

% Load the HF data.
load(welFILE,'S');

% Load the Neuquén Basin boundaries.
data=load('../data/shapes_vaca_muerta/Neuquen.txt');
Blat=data(:,2); Blon=data(:,1);

% Load in fault data.
c1=shaperead('../data/digitalizados/c1.shp');
c2=shaperead('../data/digitalizados/c2.shp');
s_norm=shaperead('../data/digitalizados/sil-norm.shp');
s_inv=shaperead('../data/digitalizados/silv_inv.shp');
s_trans=shaperead('../data/digitalizados/silv_trans.shp');
zone=shaperead('../data/digitalizados/zonas_transferencia.shp');




%%% PLOT.
GREY=[0.85,0.85,0.85];

% Make the magnitude size axis.
Rs=getMscale(EQmag)/10;

% Get plot aspect ratio.
Va=(max(latB)-min(latB))/(max(lonB)-min(lonB))/cosd(mean(latB));

% Map.
figure(3); clf; 
plot(Blon,Blat,'-r'); hold on;
plot(-68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
plot([S.Slon],[S.Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(EQmag>=Mc),EQlat(EQmag>=Mc),Rs(EQmag>=Mc),'ok','MarkerFaceColor','r');
for i=1:length(c1); plot(c1(i).X,c1(i).Y,'-k'); end
for i=1:length(c2); plot(c2(i).X,c2(i).Y,'-k'); end
for i=1:length(s_norm); plot(s_norm(i).X,s_norm(i).Y,'-k'); end
for i=1:length(s_inv); plot(s_inv(i).X,s_inv(i).Y,'-k'); end
for i=1:length(s_trans); plot(s_trans(i).X,s_trans(i).Y,'-k'); end
for i=1:length(zone); plot(zone(i).X,zone(i).Y,'-k'); end
xlabel('Longitude'); ylabel('Latitude');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
pbaspect([1 Va 1]);











%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

