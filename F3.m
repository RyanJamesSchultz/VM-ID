% Script to plot the (visual) definition of the SAF.
% Used to make Figure 3.
clear;

% Predefine some stuff.
Rmax=25;
Tmax=60;
dR= [7 15]; Pr=0.95;
dTs=[0 30]; Pt=1.00;
n=1e3;

% Make a fake operation.
OP(1).Slat=0;
OP(1).Slon=0;
OP(1).T(1)=datetime(2000,01,16,00,00,00);
OP(1).T(2)=datetime(2000,01,31,00,00,00);
OP(1).IDe=1;

% Make a fictional set of 'earthquakes'.
f=6371*pi()/180;
latEQ=linspace(0,Rmax,n)/f;
lonEQ=zeros(size(latEQ));
Rnull=lonEQ;
tEQ=  linspace(OP(1).T(1),OP(1).T(2)+Tmax,n);
Tnull=linspace(OP(1).T(1),OP(1).T(1)     ,n);

% Do the SAF processing for spatial-scores.
[S1r]=SAF(OP,latEQ,lonEQ,Tnull,    dR,  Pr,    dTs, Pt,'linear','unique'); Lr1=Rnull; Lr1(S1r.Ieq)=[S1r.Lr]; % This study.
[S2r]=SAF(OP,latEQ,lonEQ,Tnull,[5 10],  Pr,    dTs, Pt,'linear','unique'); Lr2=Rnull; Lr2(S2r.Ieq)=[S2r.Lr]; % This study (restrictive spatial-score).
[S3r]=SAF(OP,latEQ,lonEQ,Tnull,[2  1],1.00,[5 3.5],1.0,  'exp2','unique'); Lr3=Rnull; Lr3(S3r.Ieq)=[S3r.Lr]; % Lomax & Savvaidis, 2019.
[S4r]=SAF(OP,latEQ,lonEQ,Tnull,[4  2],1.00,[5  10],1.0,  'exp2','unique'); Lr4=Rnull; Lr4(S4r.Ieq)=[S4r.Lr]; % Savvaidis et al., 2020.
[S5r]=SAF(OP,latEQ,lonEQ,Tnull,[3  -1.2137],3.794,[5  -0.79063],3.5696, 'power','unique'); Lr5=Rnull; Lr5(S5r.Ieq)=[S5r.Lr]; % Ghofrani & Atkinson, 2020.
[S6r]=SAF(OP,latEQ,lonEQ,Tnull,[3  -1.2137],3.794,[1  -0.76862],1.0000, 'power','unique'); Lr6=Rnull; Lr6(S6r.Ieq)=[S6r.Lr]; % Ghofrani & Atkinson, 2021.

% Do the SAF processing for temporal-scores.
[S1t]=SAF(OP,Rnull,lonEQ,tEQ,      dR,  Pr,    dTs, Pt,'linear','unique'); Lt1=Rnull; Lt1(S1t.Ieq)=[S1t.Lt]; % This study.
[S2t]=SAF(OP,Rnull,lonEQ,tEQ,  [5 10],  Pr,    dTs, Pt,'linear','unique'); Lt2=Rnull; Lt2(S2t.Ieq)=[S2t.Lt]; % This study (restrictive spatial-score).
[S3t]=SAF(OP,Rnull,lonEQ,tEQ,  [2  1],1.00,[5 3.5],1.0,  'exp2','unique'); Lt3=Rnull; Lt3(S3t.Ieq)=[S3t.Lt]; % Lomax & Savvaidis, 2019.
[S4t]=SAF(OP,Rnull,lonEQ,tEQ,  [4  2],1.00,[5  10],1.0,  'exp2','unique'); Lt4=Rnull; Lt4(S4t.Ieq)=[S4t.Lt]; % Savvaidis et al., 2020.
[S5t]=SAF(OP,Rnull,lonEQ,tEQ,  [3  -1.2137],3.794,[5  -0.79063],3.5696, 'power','unique'); Lt5=Rnull; Lt5(S5t.Ieq)=[S5t.Lt]; % Ghofrani & Atkinson, 2020.
[S6t]=SAF(OP,Rnull,lonEQ,tEQ,  [3  -1.2137],3.794,[1  -0.76862],1.0000, 'power','unique'); Lt6=Rnull; Lt6(S6t.Ieq)=[S6t.Lt]; % Ghofrani & Atkinson, 2021.

% Make plotting axes.
T=daysdif(OP(1).T(2),tEQ);
R=latEQ*f;

% Plot Figure 3.
figure(3); clf;
% L(r).
subplot(121); hold on;
plot(R,Lr1,'-', 'Color','#0072BD','DisplayName','This Study');
plot(R,Lr2,'--','Color','#4DBEEE','DisplayName','This Study (restrictive spatial-score)');
plot(R,Lr3,'-', 'Color','#D95319','DisplayName','Lomax & Savvaidis, 2019');
plot(R,Lr4,'-', 'Color','#A2142F','DisplayName','Savvaidis et al., 2020');
plot(R,Lr5,'-', 'Color','#00FF00','DisplayName','Ghofrani & Atkinson, 2020; 2021');
%plot(R,Lr5,'DisplayName','Ghofrani & Atkinson, 2021');
plot(dR(1)*[1 1], [0 1.05],'--k','HandleVisibility','off');
xlabel('Event-well Distance (km)'); ylabel('Event-well Spatial-score');
title('L(r)');
legend();
xlim([0 max(Rmax)]); ylim([0 1.05]);
% L(t).
subplot(122); hold on;
plot(T,Lt1,'-', 'Color','#0072BD','DisplayName','This Study');
plot(T,Lt3,'-', 'Color','#D95319','DisplayName','Lomax & Savvaidis, 2019');
plot(T,Lt4,'-', 'Color','#A2142F','DisplayName','Savvaidis et al., 2020');
plot(T,Lt5,'-', 'Color','#00FF00','DisplayName','Ghofrani & Atkinson, 2020');
plot(T,Lt6,'-', 'Color','#77AC30','DisplayName','Ghofrani & Atkinson, 2021');
plot(daysdif(OP(1).T(2),[OP(1).T(1), OP(1).T(1)]), [0 1.05],'--k','HandleVisibility','off');
plot(daysdif(OP(1).T(2),[OP(1).T(2), OP(1).T(2)]), [0 1.05],'--k','HandleVisibility','off');
xlabel('Time after well shut-in (days)'); ylabel('Event-well Temporal-score');
xlim(daysdif(OP(1).T(2),[min(tEQ)-5 max(tEQ)])); ylim([0 1.05]);
title('L(t)');
legend();


