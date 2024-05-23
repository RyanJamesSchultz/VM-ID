% Script used to boostrap the spatiotemporal association filtering on the 
% Neuquén Basin data.  Used to inform Section 2.2 and make Figures 3 & S4.
clear;

% Predefine some stuff.
dM=0.3;
Mc=2.5;
dR=[7.0 15.0];
dTs=30;
Scut=1;
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/20240314_VM_quakes_all_soruces.csv';
welFILE='/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat';
Nb=5e2;

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
B(1).S=[L.S];

% Bootstrap trials.
for i=1:Nb
    i

    % Reshuffle the EQ catalogue in time.
    EQtime_i=EQtime(randperm(length(EQtime)));
    %Sj=S; [Sj.]
    
    % Do SAF and save data.
    Li=SAF(S,EQlat,EQlon,EQtime_i,dR,dTs,'unique');
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




%%% PLOT.
GREY=[0.85,0.85,0.85];

% Bootstrapping histograms.
figure(1); clf;

% Basin-scale well association percentages.
subplot(121); hold on;
histogram(100*Pi,round(0.75*sqrt(length(Pi))), 'FaceColor', GREY);
plot(100*[P1 P1],ylim(),'-r');
xlabel('Wells Associated (%)');  ylabel('Counts');

% SAF-score histograms.
subplot(122); hold on;
Sa=[B(2:end).S];
histogram(S1(S1>0),0:11,'Normalization','pdf', 'FaceColor', 'r');
histogram(Sa(Sa>0),0:11,'Normalization','pdf', 'FaceColor', GREY);
xlabel('SAF-score'); ylabel('Probability Density');
xlim([0 10]);

% Bootstrapped p-values.
figure(2); clf; hold on;
p=log10([B(2:end).p]);
histogram(p,round(2*sqrt(length(p))));
plot(log10(0.05)*[1 1],ylim(),'--k');
xlabel('log_{10}(p-value)'); ylabel('Counts');
xlim([min(p)-0.25 0]);

% Report the SAF-score KS-test results.
1-10.^mean(p)







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