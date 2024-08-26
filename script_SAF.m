% Script used to do the spatiotemporal association filtering on the Neuquén Basin data.  
% Used to inform Section 3.1 and make Figures 4 & S6-S10.
clear;

% Predefine some stuff.
Mc=2.5; dM=0.3;
dR= [7 15]; Pr=0.95;
dTs=[0 30]; Pt=1.00;
decay_flag='linear'; assoc_type='unique';
Scut1=1;
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/CAT.mat';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';

% Boundary details and for filtering.
latB=-1;
lonB=-1;
depB=-1;
%tB=[datetime(2019,03,01) datetime(2019,04,01)];
tB=-1;
mB=-1;

% Load the EQ catalogue data.
load(catFILE,'EQlat','EQlon','EQdep','EQtime','EQmag');
[EQlat,EQlon,EQdep,EQmag,EQtime]=filtEQ(EQlat,EQlon,EQdep,EQmag,EQtime,  latB,lonB,tB,mB);

% Load the HF data.
load(welFILE,'S');
S=filtHF(S,latB,lonB,tB);

% Load the Neuquén Basin boundaries.
data=load('/Users/rschultz/Desktop/VacaMuerta/data/shapes_vaca_muerta/Neuquen.txt');
Blat=data(:,2); Blon=data(:,1);

% GR-MFD stats.
[b,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(EQmag, Mc,dM);
po=[-b,a];
Mgr_fit=[Mc, max(EQmag)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Associate HF-EQ datasets together.
L=SAF(S,EQlat,EQlon,EQtime,dR,Pr,dTs,Pt,decay_flag,assoc_type); Scut2=5.5; % This study.
%L=SAF(S,EQlat,EQlon,EQtime,[5 10],Pr,dTs,Pt,decay_flag,assoc_type); Scut2=3.7; % This study (restrictive spatial-score).
%L=SAF(S,EQlat,EQlon,EQtime,[2 1],0.99,[5 3.5],0.99,'exp2','unique'); Scut2=2.0; % Lomax & Savvaidis, 2019.
%L=SAF(S,EQlat,EQlon,EQtime,[4 2],0.99,[5  10],0.99,'exp2','unique'); Scut2=3.4; % Savvaidis et al., 2020.
%L=SAF(S,EQlat,EQlon,EQtime,[3 -1.2137],3.794,[5 -0.79063],3.5696,'power','unique'); Scut2=3.8; % Ghofrani & Atkinson, 2020.
%L=SAF(S,EQlat,EQlon,EQtime,[3 -1.2137],3.794,[1 -0.76862],1.0000,'power','unique'); Scut2=3.8; % Ghofrani & Atkinson, 2021.

% Find the wells and EQs that pass cut-off criteria.
Iop1=find([L.S]>=Scut1);
Iop2=find([L.S]>=Scut2);
Ieq1=unique(vertcat(L(Iop1).Ieq));
Ieq2=unique(vertcat(L(Iop2).Ieq));

% Find the HF wells with the highest SAF-scores.
[~,Itw]=sort([L.S],'descend');

% Report the percent association results.
disp('% Asscociations')
(length(Iop1)/length(L))*100
(length(Ieq1)/length(EQtime))*100

% Report the top SAF-scoring wells.
[~,Is]=sort([L(Iop2).S],'descend');
S([L(Iop2(Is)).Ihf]).UWI





%%% PLOT.
GREY=[0.85,0.85,0.85];

% Make the magnitude size axis.
Rs=getMscale(EQmag)/10;

% MvT
figure(1); clf; 
ax1=subplot(211); hold on;
for i=1:length(S)
    plot([S(i).T],log10(S(i).Vt)*[1 1],'-db','MarkerFaceColor','c');
end
scatter(EQtime,EQmag,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQtime(EQmag>=Mc),EQmag(EQmag>=Mc),Rs(EQmag>=Mc),'ok','MarkerFaceColor','r');
xlabel('Time'); ylabel('Magntidue');
ax2=subplot(212); hold on;
for i=1:length(S)
    if(L(i).S>=Scut1)
        plot([S(i).T],log10(S(i).Vt)*[1 1],'-db','MarkerFaceColor','c');
    else
        %plot([S(i).T],log10(S(i).Vt)*[1 1],'-db');
    end
end
scatter(EQtime,EQmag,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQtime(Ieq1),EQmag(Ieq1),Rs(Ieq1),'ok','MarkerFaceColor','r');
xlabel('Time'); ylabel('Magntidue');
linkaxes([ax1 ax2],'x');

% The GR-FMD.
figure(2); clf;
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
plot(Mc*[1 1],ylim,'--k');
xlabel('Magnitude (M_L)'); ylabel('Count');
xlim([min(Mgr)-dM/2 max(Mgr)+dM/2]); ylim([0.7 1.3*max(Ngr)]);

% Map.
figure(3); clf; 
ax1=subplot(121); hold on;
plot(Blon,Blat,'-k');
plot(-68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
plot([S.Slon],[S.Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(EQmag>=Mc),EQlat(EQmag>=Mc),Rs(EQmag>=Mc),'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
ax2=subplot(122); hold on;
plot(Blon,Blat,'-k');
plot( -68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
plot([S.Slon],[S.Slat],'bd');
plot([S(Iop1).Slon],[S(Iop1).Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(Ieq1),EQlat(Ieq1),Rs(Ieq1),'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
linkaxes([ax1 ax2],'xy');

% Figure 4.
latB=[-39.2 -37.3]; lonB=[-69.75 -67.55];
Va=(max(latB)-min(latB))/(max(lonB)-min(lonB))/cosd(mean(latB));
figure(4); clf;
% Full Map
subplot(4,5,[1,2,6,7]); hold on;
plot(Blon,Blat,'-k');
plot(-68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
plot([S.Slon],[S.Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(EQmag>=Mc),EQlat(EQmag>=Mc),Rs(EQmag>=Mc),'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
pbaspect([1 Va 1]);
% Full MvT.
subplot(4,5,[3,4,5,8,9,10]); hold on;
for i=1:length(S)
    plot([S(i).T],log10(S(i).Vt)*[1 1],'-db','MarkerFaceColor','c');
end
scatter(EQtime,EQmag,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQtime(EQmag>=Mc),EQmag(EQmag>=Mc),Rs(EQmag>=Mc),'ok','MarkerFaceColor','r');
xlabel('Time'); ylabel('Magntidue');
xlim([min([S.T]) max([S.T])]); ylim([0 6]);
% SAF Map.
subplot(4,5,[11,12,16,17]); hold on;
plot(Blon,Blat,'-k');
plot( -68.06,-38.96,'ks','MarkerFaceColor','k');
plot([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],'bs','MarkerFaceColor','b');
plot([S.Slon],[S.Slat],'bd');
plot([S(Iop1).Slon],[S(Iop1).Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(Ieq1),EQlat(Ieq1),Rs(Ieq1),'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
pbaspect([1 Va 1]);
% SAF MvT.
subplot(4,5,[13,14,15,18,19,20]); hold on;
for i=1:length(S)
    if(L(i).S>=Scut1)
        plot([S(i).T],log10(S(i).Vt)*[1 1],'-db','MarkerFaceColor','c');
    else
        %plot([S(i).T],log10(S(i).Vt)*[1 1],'-db');
    end
end
scatter(EQtime,EQmag,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQtime(Ieq1),EQmag(Ieq1),Rs(Ieq1),'ok','MarkerFaceColor','r');
xlabel('Time'); ylabel('Magntidue');
xlim([min([S.T]) max([S.T])]); ylim([0 6]);

% ID finder.
figure(5); clf; 
ax3=gca;
plot3(Blon,Blat,zeros(size(Blon)),'-k'); hold on;
plot3(-68.06,-38.96,length(S),'ks','MarkerFaceColor','k'); 
plot3([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],length(S)*[1 1 1],'bs','MarkerFaceColor','b');
scatter3([S.Slon],[S.Slat],[S.IDe],'bd','MarkerFaceColor','c');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Well ID');
linkaxes([ax1 ax2 ax3],'xy');





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