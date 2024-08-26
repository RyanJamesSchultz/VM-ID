% Script to make Figure 2 - the Gutenberg-Richter magnitude-frequency distribution.
clear;

% Predefine some stuff.
dM=0.3;
Mc=2.5;
GREY=[0.85,0.85,0.85];
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/CAT.mat';

% Boundary polygons and stuff.
latB=-1;
lonB=-1;
depB=-1;
%tB=[datetime(2017,01,01) datetime(2018,01,01)];
tB=-1;
mB=-1;

% Load the EQ catalogue data.
load(catFILE,'EQlat','EQlon','EQdep','EQtime','EQmag');
[EQlat,EQlon,EQdep,EQmag,EQtime]=filtEQ(EQlat,EQlon,EQdep,EQmag,EQtime,  latB,lonB,tB,mB);

% GR-MFD stats.
[b,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(EQmag, Mc,dM);
po=[-b,a];
Mgr_fit=[Mc, max(EQmag)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Plot the GR-FMD.
figure(2); clf;
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
xlim([min(Mgr)-dM/2 max(Mgr)+dM/2]); ylim([0.7 1.3*max(Ngr)]);
plot(Mc*[1 1],ylim,'--k');
xlabel('Magnitude (M)'); ylabel('Count');




%%%% SUBROUNTINES.

% Spatiotemporally filter the EQ catalogue.
function [lat,lon,dep,M,T]=filtEQ(lat,lon,dep,M,T,  lat_L,lon_L,T_L,M_L)
  
  % Filter spatially (lateral).
  if(lat_L~=-1)
      I=inpolygon(lon,lat,lon_L,lat_L);
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
  end
  
  % Filter temporally.
  if(length(T_L)==2)
      I=(T>=min(T_L))&(T<=max(T_L));
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
  end
  
  % Filter by magnitudes.
  if(M_L~=-1)
      I=(M>=min(M_L))&(M<=max(M_L));
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
  end
  
end
