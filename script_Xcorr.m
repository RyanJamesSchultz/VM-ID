% Script that will perform the cross-correlation reshuffling analysis on the top SAF-scoring events.  
% Used to inform Section 3.2 and make Figures 6 & S26-SS31.
clear;

% Predefine some stuff.
Mc=2.5; dM=0.3;
dR= [7 15]; Pr=0.95;
dTs=[0 30]; Pt=1.00;
decay_flag='linear'; assoc_type='unique';
dTx=1; Nx=1e5; conf=[0.6827 0.9545 0.9973 0.9999];
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/CAT.mat';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';

% UWI list for filtering of cases.
UWIs={'YPF.Nq.LACh-45(h)','TPT.Nq.FP-1316(h)','TPT.Nq.FP-1427(h)','YPF.Nq.LACh-385(h)','YPF.Nq.LACh-161(h)','TAU.Nq.APIg-112(h)','SHE.Nq.BAñ-1015(h)','TPT.Nq.FP-1195(h)','TPT.Nq.FP-1425(h)','YPF.Nq.AdCh-1186(h)','SHE.Nq.BAñ-1011(h)'};
UWIs=UWIs(11);

% Boundary details and for filtering.
latB=-1;
lonB=-1;
depB=-1;
tB=-1;
mB=-1;

% Load the EQ catalogue data.
load(catFILE,'EQlat','EQlon','EQdep','EQtime','EQmag');
[EQlat,EQlon,EQdep,EQmag,EQtime]=filtEQ(EQlat,EQlon,EQdep,EQmag,EQtime,  latB,lonB,tB,mB);

% Get the HF data.
load(welFILE,'S');
S=filtHF(S,latB,lonB,tB);

% Associate HF-EQ datasets together.
L=SAF(S,EQlat,EQlon,EQtime,dR,Pr,dTs,Pt,decay_flag,assoc_type);

% Filter the HF data to just the subset of wells.
Iop=find(cellfun(@(x) ismember(x, UWIs), {S.UWI}));
S=S(Iop);
L=L(Iop);

% Filter the EQ data to just the subset of associated events.
Ieq=unique(vertcat(L.Ieq));
EQlat=EQlat(Ieq);
EQlon=EQlon(Ieq);
EQdep=EQdep(Ieq);
EQmag=EQmag(Ieq);
EQtime=EQtime(Ieq);

% Make time axis (volume).
Ts=min([S.T])-days(dTs(2));
Te=max([S.T])+days(dTs(2));
Tt=Ts:days(dTx):Te;
Tc=daysdif(Ts,Tt);

% Make the volume axis.
Vt=zeros(size(Tc));
for i=1:length(S)
    r=S(i).V(2)/(daysdif(S(i).T(1),S(i).T(2))*24*60);
    Vt((Tt>=S(i).T(1))&(Tt<S(i).T(2)))=Vt((Tt>=S(i).T(1))&(Tt<S(i).T(2)))+r;
end

% Make EQ rate & time axes.
Teq=daysdif(Ts,EQtime);
Neq1=histcounts(Teq(EQmag>=Mc),[Tc-dTx/2, Tc(end)+dTx/2]);
Neq2=histcounts(Teq,           [Tc-dTx/2, Tc(end)+dTx/2]);

% Do the reshuffling cross-correlation analysis.
[Tx,CC,CC_conf]=ReshuffleCorr(Tc, Neq1, Vt, Nx, conf);





%%% Plot.
GREY=[0.85,0.85,0.85];
PURP=[133,57,227]/256;

% Timeseries and xcorr plots.
figure(6); clf;
subplot(121);
plot(Tt,Vt,'-b'); hold on;
bar(Tt,Neq2, 'FaceColor', GREY);
bar(Tt,Neq1, 'FaceColor', 'r');
xlabel('Time'); ylabel('Earthquake Counts');
if length(UWIs)==1, title(UWIs), end
subplot(122);
plot(Tx,CC,'-','Color',PURP); hold on;
plot(Tx, CC_conf(:,1),'-','Color','k');
plot(Tx, CC_conf(:,2),'-','Color','k');
plot(Tx, CC_conf(:,3),'-','Color','k');
plot(Tx, CC_conf(:,4),'-','Color','k');
xlabel('Lag Time (days)'); ylabel('Cross-correlation Amplitude');
xlim([-30 +40]);





%%%% SUBROUNTINES.

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