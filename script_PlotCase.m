% Script to plot EQ and HF data around an SAF identified case.
% Used to inform Section 3.1 and make Figures S16-S25.
clear;

% Predefine some stuff.
Mc=2.5;
dX=0.15;
dR= [7 15]; Pr=0.95;
dTs=[0 30]; Pt=1.00;
decay_flag='linear'; assoc_type='unique';
Scut1=1;
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/CAT.mat';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';

% UWI list for filtering of cases.
UWIs={'YPF.Nq.LACh-45(h)','TPT.Nq.FP-1316(h)','TPT.Nq.FP-1427(h)','YPF.Nq.LACh-385(h)','YPF.Nq.LACh-161(h)','TAU.Nq.APIg-112(h)','SHE.Nq.BAñ-1015(h)','TPT.Nq.FP-1195(h)','TPT.Nq.FP-1425(h)','YPF.Nq.AdCh-1186(h)','SHE.Nq.BAñ-1011(h)'};
UWIs=UWIs(11);

% Load the EQ & HF data.
load(catFILE,'EQlat','EQlon','EQdep','EQtime','EQmag');
load(welFILE,'S');

% Load the Neuquén Basin boundaries.
data=load('/Users/rschultz/Desktop/VacaMuerta/data/shapes_vaca_muerta/Neuquen.txt');
Blat=data(:,2); Blon=data(:,1);

% Filter the HF data to just the subset of well UWIs.
Iop=find(cellfun(@(x) ismember(x, UWIs), {S.UWI}));
Sf=S(Iop);

% Boundary details and for filtering.
xm=min(arrayfun(@(s) s.Slat,Sf)); xM=max(arrayfun(@(s) s.Slat,Sf));
ym=min(arrayfun(@(s) s.Slon,Sf)); yM=max(arrayfun(@(s) s.Slon,Sf));
latB=[xm-dX xM+dX xM+dX xm-dX];
lonB=[yM+dX yM+dX ym-dX ym-dX];
depB=-1;
tB=[min(arrayfun(@(s) s.T(1),Sf))-10 max(arrayfun(@(s) s.T(2),Sf))+60];
mB=-1;

% Filter the datasets.
[EQlat,EQlon,EQdep,EQmag,EQtime]=filtEQ(EQlat,EQlon,EQdep,EQmag,EQtime,  latB,lonB,tB,mB);
S=filtHF(S,latB,lonB,tB);

% Associate HF-EQ datasets together.
L=SAF(Sf,EQlat,EQlon,EQtime,dR,Pr,dTs,Pt,decay_flag,assoc_type); Scut2=5.5; % This study.
%L=SAF(S,EQlat,EQlon,EQtime,[5 10],Pr,dTs,Pt,decay_flag,assoc_type); Scut2=3.7; % This study (restrictive spatial-score).
%L=SAF(S,EQlat,EQlon,EQtime,[2 1],0.99,[5 3.5],0.99,'exp2','unique'); Scut2=2.0; % Lomax & Savvaidis, 2019.
%L=SAF(S,EQlat,EQlon,EQtime,[4 2],0.99,[5  10],0.99,'exp2','unique'); Scut2=3.4; % Savvaidis et al., 2020.
%L=SAF(S,EQlat,EQlon,EQtime,[3 -1.2137],3.794,[5 -0.79063],3.5696,'power','unique'); Scut2=3.8; % Ghofrani & Atkinson, 2020.
%L=SAF(S,EQlat,EQlon,EQtime,[3 -1.2137],3.794,[1 -0.76862],1.0000,'power','unique'); Scut2=3.8; % Ghofrani & Atkinson, 2021.

% Find the wells and EQs that pass cut-off criteria.
Iop1=find([L.S]>=Scut1);
Ieq1=unique(vertcat(L(Iop1).Ieq));

% Compute the cluster centroid.
latC=mean(EQlat(Ieq1));
lonC=mean(EQlon(Ieq1));

% Compute distances from the well location and cluster centroid.
dRa=Geoid_Distance(Sf(1).Slat,Sf(1).Slon, EQlat(Ieq1),EQlon(Ieq1), 'elliptical')*6371*pi()/180; % km.
dRp=Geoid_Distance(      latC,      lonC, EQlat(Ieq1),EQlon(Ieq1), 'elliptical')*6371*pi()/180; % km.

% Report the cluster's spatial accuracy and precision.
std(dRa)
std(dRp)





%%% PLOT.
GREY=[0.85,0.85,0.85];

% Make the magnitude size axis.
Rs=getMscale(EQmag)/5;

% Figure S15-S26.
%latB=[-39.2 -37.3]; lonB=[-69.75 -67.55];
Va=(max(latB)-min(latB))/(max(lonB)-min(lonB))/cosd(mean(latB));
figure(516); clf;
% SAF Map.
subplot(1,5,[1 2]); hold on;
plot(Blon,Blat,'-k');
plot( -68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
plot([S.Slon],[S.Slat],'bd');
plot([Sf.Slon],[Sf.Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(Ieq1),EQlat(Ieq1),Rs(Ieq1),'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
pbaspect([1 Va 1]);
% SAF MvT.
subplot(1,5,[3 4 5]); hold on;
for i=1:length(S)
    plot([S(i).T],log10(S(i).Vt)*[1 1],'-db');
end
for i=1:length(Sf)
    plot([Sf(i).T],log10(Sf(i).Vt)*[1 1],'-db','MarkerFaceColor','c');
end
scatter(EQtime,EQmag,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQtime(Ieq1),EQmag(Ieq1),Rs(Ieq1),'ok','MarkerFaceColor','r');
xlabel('Time'); ylabel('Magntidue');
xlim([tB(1) tB(2)]); ylim([0 6]);
if length(UWIs)==1, title(UWIs), end





%%%% SUBROUNTINES.

% Get the size of rupture from event magnitude.
function [Rs]=getMscale(ML)
  %Mw=0.754*ML+0.884;% ML to Mw [Yenier, 2017; Ross et al., 2016].
  Mw=ML;
  Mo=10.^(1.5*Mw+9.1); % Mw to Mo (Nm).
  Rs=nthroot((7/16)*(Mo/3e6),3); % Mo to radius (m), assuming a stress dop of 3 MPa.
end

% Spatiotemporally filter the HF dataset.
function [S]=filtHF(S,lat_L,lon_L,T_L)
  
  % Filter spatially (lateral).
  if(lat_L~=-1)
      I=inpolygon([S.Slon],[S.Slat],lon_L,lat_L);
      S=S(I);
  end
  
  % Filter temporally.
  if(length(T_L)==2)
      Ts=arrayfun(@(S) S.T(1),S);
      Te=arrayfun(@(S) S.T(2),S);
      I=(Te>=min(T_L))&(Ts<=max(T_L));
      S=S(I);
  end
  
end

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

