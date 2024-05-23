% Script to make Figure S2 - the Gutenberg-Richter magnitude-frequency distribution.
clear;

% Predefine some stuff.
dM=0.3;
Mc=2.6;
GREY=[0.85,0.85,0.85];

% Boundary polygons and stuff.
latB=-1;
lonB=-1;
depB=-1;
%tB=[datetime(2017,01,01) datetime(2018,01,01)];
tB=-1;
mB=-1;

% Load in catalogue data.
catFILE='/Users/rschultz/Desktop/VacaMuerta/data/20240314_VM_quakes_all_soruces.csv';
[EQlat,EQlon,EQdep,EQtime,EQmag]=parseINPRES(catFILE,latB,lonB,depB,tB,mB);

% GR-MFD stats.
[b,b_err,a,R2,~,Mgr,Ngr,ngr]=Bval(EQmag, Mc,dM);
po=[-b,a];
Mgr_fit=[Mc, max(EQmag)];
Ngr_fit=10.^polyval(po,Mgr_fit);

% Plot the GR-FMD.
figure(52); clf;
semilogy(Mgr, Ngr, 'o', 'Color', 'k'); hold on;
bar(Mgr,ngr, 'FaceColor', GREY);
semilogy(Mgr_fit, Ngr_fit, '-', 'Color', 'black');
xlim([min(Mgr)-dM/2 max(Mgr)+dM/2]); ylim([0.7 1.3*max(Ngr)]);
plot(Mc*[1 1],ylim,'--k');
xlabel('Magnitude (M)'); ylabel('Count');
