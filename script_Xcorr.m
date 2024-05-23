% Script that will perform the cross-correlation reshuffling analysis on
% a list of wells and SAF-associated events.  Used to inform Section 2.3
% and make Figures 4 & S5-S7.
clear;

% Predefine some stuff.
dM=0.3;
Mc=2.5;
dR=[7.0 15.0];
dTs=30;
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/20240314_VM_quakes_all_soruces.csv';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';
UWIs={'TPT.NQ.FP-1427(H)','TPT.NQ.FP-1316(H)','SHE.NQ.BAN-1015(H)','YPF.NQ.LACH-161(H)','TPT.NQ.FP-1425(H)','YPF.NQ.LACH-157(H)','SHE.NQ.BAN-1011(H)','VIS.NQ.BPO-2394(H)','TPT.NQ.FP-1195(H)','YPF.NQ.RDMN-229(H)','YPF.NQ.RDMN-48(H)','PLU.NQ.LCA-3086(H)'};
%UWIs=UWIs(1);
dTx=1;
Nx=1e5;
conf=[0.6827 0.9545 0.9973 0.9999];

% Boundary polygons and stuff.
latB=-1;
lonB=-1;
depB=-1;
tB=-1;
mB=-1;

% Get the EQ catalogue data.
[EQlat,EQlon,EQdep,EQtime,EQmag]=parseINPRES(catFILE,latB,lonB,depB,tB,mB);

% Get the HF data.
load(welFILE,'S');
S=filtHF(S,latB,lonB,tB);

% Associate HF-EQ datasets together.
L=SAF(S,EQlat,EQlon,EQtime,dR,dTs,'unique');


% Filter the HF data to just the subset of wells.
Iop=find(cellfun(@(x) ismember(x, UWIs), {S.UWI}));
%Iop=find([L.S]>=Seq);
%%%SORT TO BE IN THE SAME ORDER AS THE INPUT?
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
Ts=min([S.T])-days(dTs);
Te=max([S.T])+days(dTs);
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
figure(4); clf;
subplot(121);
plot(Tt,Vt,'-b'); hold on;
bar(Tt,Neq2, 'FaceColor', GREY);
bar(Tt,Neq1, 'FaceColor', 'r');
xlabel('Time'); ylabel('Earthquake Counts');
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