% Script used to boostrap the spatiotemporal association filtering on the Neuquén Basin data.  
% Used to inform Section 3.1 and make Figures 5, S5, & S11-S15.
clear;

% Predefine some stuff.
Mc=2.5; dM=0.3;
dR= [7 15]; Pr=0.95;
dTs=[0 30]; Pt=1.00;
decay_flag='linear'; assoc_type='unique';
Scut=1;
Nb=10;
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/CAT.mat';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';

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

% Associate HF-EQ datasets together via SAF.
L=SAF(S,EQlat,EQlon,EQtime,dR,Pr,dTs,Pt,decay_flag,assoc_type); % This study.
%L=SAF(S,EQlat,EQlon,EQtime,[2 1],0.99,[5 3.5],0.99,'exp2','unique'); % Lomax & Savvaidis, 2019.
%L=SAF(S,EQlat,EQlon,EQtime,[4 2],0.99,[5  10],0.99,'exp2','unique'); % Savvaidis et al., 2020.
%L=SAF(S,EQlat,EQlon,EQtime,[3 -1.2137],3.794,[5 -0.79063],3.5696,'power','unique'); % Ghofrani & Atkinson, 2020.
%L=SAF(S,EQlat,EQlon,EQtime,[3 -1.2137],3.794,[1 -0.76862],1.0000,'power','unique'); % Ghofrani & Atkinson, 2021.
B(1).S=[L.S];

% Bootstrap trials.
for i=1:Nb
    i
    
    % Reshuffle the EQ catalogue in time.
    EQtime_i=EQtime(randperm(length(EQtime)));
    
    % Do the SAF processing and save data.
    Li=SAF(S,EQlat,EQlon,EQtime_i,dR,Pr,dTs,Pt,decay_flag,assoc_type); % This study.
    %Li=SAF(S,EQlat,EQlon,EQtime_i,[2 1],0.99,[5 3.5],0.99,'exp2','unique'); % Lomax & Savvaidis, 2019.
    %Li=SAF(S,EQlat,EQlon,EQtime_i,[4 2],0.99,[5  10],0.99,'exp2','unique'); % Savvaidis et al., 2020.
    %Li=SAF(S,EQlat,EQlon,EQtime_i,[3 -1.2137],3.794,[5 -0.79063],3.5696,'power','unique'); % Ghofrani & Atkinson, 2020.
    %Li=SAF(S,EQlat,EQlon,EQtime_i,[3 -1.2137],3.794,[1 -0.76862],1.0000,'power','unique'); % Ghofrani & Atkinson, 2021.
    B(i+1).S=[Li.S];
end

% Compute the percentage of association.
S1=[B(1).S]; P1=sum(S1>Scut)/length(S1);

% Loop over all of the bootstrap trials.
for i=2:length(B)
    
    % Compute the (reshuffled) percentage of association.
    Si=[B(i).S]; Pi(i-1)=sum(Si>Scut)/length(Si);
    
    % Do KS-test on (non-zero) SAF-scores.
    [~,p]=kstest2(S1(S1>0),Si(Si>0));
    B(i).p=(p);
end

% Save the workspace.
save('SAFbs.mat');





%%
%clear; load('SAFbs_Sc24.mat');

%%% PLOT.
GREY=[0.85,0.85,0.85];

% Bootstrapped histograms.
figure(5); clf;

% Basin-scale well association percentages.
subplot(121); hold on;
histogram(100*Pi,round(0.75*sqrt(length(Pi))), 'FaceColor', GREY);
plot(100*[P1 P1],ylim(),'-r');
xlabel('Wells Associated (%)');  ylabel('Counts');

% SAF-score histograms.
subplot(122); hold on;
Sa=[B(2:end).S];
histogram(S1(S1>0),0:15,'Normalization','pdf', 'FaceColor', 'r');
histogram(Sa(Sa>0),0:15,'Normalization','pdf', 'FaceColor', GREY);
xlabel('SAF-score'); ylabel('Probability Density');
set(gca, 'YScale', 'log');
xlim([0 14]);

% Bootstrapped p-values.
figure(54); clf; hold on;
p=log10([B(2:end).p]);
histogram(p,round(2*sqrt(length(p))));
plot(log10(0.05)*[1 1],ylim(),'--k');
xlabel('log_{10}(p-value)'); ylabel('Counts');
xlim([min(p)-0.25 0]);

% Report the percent association results.
disp('% Asscociations')
P1*100
[mean(Pi) median(Pi) std(Pi) max(Pi)]*100

% Report the SAF-score & false-positive results.
disp('SAF-scores')
max(S1)
[mean(arrayfun(@(s) max(s.S) ,B(2:end))) std(arrayfun(@(s) max(s.S) ,B(2:end))) max(Sa)]

% Report the KS-test results.
disp('KS-test')
(1-10.^mean(p))*100





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