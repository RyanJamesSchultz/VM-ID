% Script used to do the spatiotemporal association filtering on the Neuquén
% Basin data.  Used to inform Section 2.2 and make Figure 2.
clear;

% Predefine some stuff.
dM=0.3;
Mc=2.5;
dR=[7.0 15.0];
dTs=30;
Scut=1;
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/20240314_VM_quakes_all_soruces.csv';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';

% Boundary polygons and stuff.
latB=-1;
lonB=-1;
depB=-1;
%tB=[datetime(2017,01,01) datetime(2018,01,01)];
tB=-1;
mB=-1;

% Load the EQ catalogue data.
[EQlat,EQlon,EQdep,EQtime,EQmag]=parseINPRES(catFILE,latB,lonB,depB,tB,mB);

% Load the HF data.
load(welFILE,'S');
S=filtHF(S,latB,lonB,tB);

% Load the Neuquén Basin boundaries.
data=load('../data/shapes_vaca_muerta/Neuquen.txt');
Blat=data(:,2); Blon=data(:,1);

% GR-MFD stats.
[b,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(EQmag, Mc,dM);
po=[-b,a];
Mgr_fit=[Mc, max(EQmag)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Associate HF-EQ datasets together.
L=SAF(S,EQlat,EQlon,EQtime,dR,dTs,'unique');
Iop=find([L.S]>=Scut);
Ieq=unique(vertcat(L(Iop).Ieq));

% Find the HF wells with the highest association scores.
[~,Ihf]=sort([L.S],'descend');

% Report the top five SAF-scoring wells.
S([L(Ihf(1:5)).Ihf]).UWI




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
    if(L(i).S>=Scut)
        plot([S(i).T],log10(S(i).Vt)*[1 1],'-db','MarkerFaceColor','c');
    else
        %plot([S(i).T],log10(S(i).Vt)*[1 1],'-db');
    end
end
scatter(EQtime,EQmag,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQtime(Ieq),EQmag(Ieq),Rs(Ieq),'ok','MarkerFaceColor','r');
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
plot([S(Iop).Slon],[S(Iop).Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(Ieq),EQlat(Ieq),Rs(Ieq),'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
linkaxes([ax1 ax2],'xy');

% ID finder.
figure(4); clf; 
ax3=gca;
plot3(Blon,Blat,zeros(size(Blon)),'-k'); hold on;
plot3(-68.06,-38.96,length(S),'ks','MarkerFaceColor','k'); 
plot3([-68.818687 -68.556231 -68.851516],[-38.548934 -38.601331 -39.371536],length(S)*[1 1 1],'bs','MarkerFaceColor','b');
scatter3([S.Slon],[S.Slat],[S.IDe],'bd','MarkerFaceColor','c');
xlabel('Longitude'); ylabel('Latitude'); zlabel('Well ID');
linkaxes([ax1 ax2 ax3],'xy');

% Figure 2.
latB=[-39.2 -37.3]; lonB=[-69.75 -67.55];
Va=(max(latB)-min(latB))/(max(lonB)-min(lonB))/cosd(mean(latB));
figure(5); clf;
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
plot([S(Iop).Slon],[S(Iop).Slat],'bd','MarkerFaceColor','c');
scatter(EQlon,EQlat,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQlon(Ieq),EQlat(Ieq),Rs(Ieq),'ok','MarkerFaceColor','r');
xlabel('Longitude'); ylabel('Latitude'); %zlabel('Depth (km)');
ylim([min(latB) max(latB)]); xlim([min(lonB) max(lonB)]);
pbaspect([1 Va 1]);
% SAF MvT.
subplot(4,5,[13,14,15,18,19,20]); hold on;
for i=1:length(S)
    if(L(i).S>=Scut)
        plot([S(i).T],log10(S(i).Vt)*[1 1],'-db','MarkerFaceColor','c');
    else
        %plot([S(i).T],log10(S(i).Vt)*[1 1],'-db');
    end
end
scatter(EQtime,EQmag,Rs,'ok','MarkerFaceColor',GREY);
scatter(EQtime(Ieq),EQmag(Ieq),Rs(Ieq),'ok','MarkerFaceColor','r');
xlabel('Time'); ylabel('Magntidue');
xlim([min([S.T]) max([S.T])]); ylim([0 6]);




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
      I=inpolygon([S.Wlon],[S.Wlat],lon_L,lat_L);
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