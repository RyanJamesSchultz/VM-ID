% Script to plot the visual definition of the SAF.
% Used to make Figure S3.
clear;

% Predefine some stuff.
dR=[7.0 15.0];
dTs=30;
f=6371*pi()/180;

% Make a fake operation.
OP(1).Slat=0;
OP(1).Slon=0;
OP(1).T(1)=datetime(2000,01,01,00,00,00);
OP(1).T(2)=datetime(2000,01,31,00,00,00);
OP(1).IDe=1;

% Make a fictional set of 'earthquakes'.
latEQ=linspace(0,max(dR),100)/f;
lonEQ=zeros(size(latEQ));
tEQ=linspace(OP(1).T(1),OP(1).T(2)+dTs,100);

% Do the SAF.
[S]=SAF(OP,latEQ,lonEQ,tEQ,dR,dTs,'unique');


% Plot.
figure(53); clf;
subplot(121); hold on;
plot(latEQ*f,[S.Ps]);
plot(dR(1)*[1 1], [0 1.05],'--k');
xlabel('Distance (km)'); ylabel('Spatial Component of Single-Event SAF-score');
xlim([0 max(dR)+1]); ylim([0 1.05]);
subplot(122); hold on;
plot(daysdif(OP(1).T(1),tEQ),[S.Pt]);
plot( daysdif(OP(1).T(1),[OP(1).T(1), OP(1).T(1)]), [0 1.05],'--k');
plot( daysdif(OP(1).T(1),[OP(1).T(2), OP(1).T(2)]), [0 1.05],'--k');
xlabel('Time (days)'); ylabel('Temporal Component of Single-Event SAF-score');
xlim(daysdif(OP(1).T(1),[min(tEQ)-1 max(tEQ)+1])); ylim([0 1.05]);