% Script that will parse a VM-HF csv file into my Matlab structure.
clear;

% Get the well header information.
Th=readtable('/Users/rschultz/Desktop/VacaMuerta/data/20240314_VM_fracture_wells.csv');

% Get the directional-survey, injection, and production data.
%Td=readtable('DirSurv.csv');
%Ti=readtable('Inj2.csv');
%Tp=readtable('Prod2.csv');
Td=[]; Ti=[]; Tp=[];

% Define the data strucutre.
S=struct('UWI',[],'IDe',[],'type',[],'Slat',[],'Slon',[],'Lh',[],'Ns',[],'T',[],'V',[],'P',[],'Vt',[],'Wlat',[],'Wlon',[],'Wtvd',[]);

% Loop over all of the well IDs.
for i =1:size(Th,1)
    
    % Stuff the header information into the injection data structure.
    uwi=Th{i,3};
    type=Th{i,17};
    S(i).IDe=i;
    S(i).UWI=uwi{1};
    S(i).type=type{1};
    S(i).Slat=Th{i,2};
    S(i).Slon=Th{i,1};
    S(i).Lh=Th{i,6};
    S(i).Ns=Th{i,7};
    
    % Get the corresponding directional survey data.
    Id=Td==S(i).IDe;
    if(~isempty(find(Id,1)))
        Wlat=Td{Id,4};
        Wlon=Td{Id,5};
        Wtvd=Td{Id,6}*0.3048/1000; % ft to km.
        Wmd=Td{Id,7};
        [~,I]=sort(Wmd);
        S(i).Wlat=Wlat(I);
        S(i).Wlon=Wlon(I);
        S(i).Wtvd=Wtvd(I);
    else
        %S(i).Wlat=[S(i).Slat; Th{i,14}];
        %S(i).Wlon=[S(i).Slon; Th{i,15}];
        %S(i).Wtvd=[0        ; Th{i,16}]*0.3048/1000; % ft to km.
    end
    
    % Get the corresponding injection data.
    Ii=Ti==S(i).IDe;
    if(~isempty(find(Ii,1)))
        T=datenum(Ti{Ii,3});
        P=Ti{Ii,4}*0.0068947573; % psi to MPa
        V=Ti{Ii,8}*0.16; % bbl to m³.
        V(isnan(V))=0;
        S(i).T=T;
        S(i).V=V;
        S(i).P=P;
        S(i).Vt=sum(V);
    else
        S(i).T=[datetime(num2str(Th{i,11}),'InputFormat','yyyyMMdd') datetime(num2str(Th{i,12}),'InputFormat','yyyyMMdd')];
        S(i).V=[0 Th{i,9}]; % m³.
        S(i).P=[Th{i,20}]*0.0068947573; % MPa
        S(i).Vt=Th{i,9}; % m³.
    end

    % Get the corresponding production data.
    Ip=Tp==S(i).IDe;
    if(~isempty(find(Ip,1)))
        tr=datenum(Tp{Ip,5});
        Ro=Tp{Ip,11}*0.16;  % bbl to m³.
        Rg=Tp{Ip,12}*28.32; % MCF to m³.
        Rw=Tp{Ip,13}*0.16;  % bbl to m³.
        Ro(isnan(Ro))=0; Rg(isnan(Rg))=0; Rw(isnan(Rw))=0;
        S(i).tr=tr;
        S(i).Ro=Ro;
        S(i).Rg=Rg;
        S(i).Rw=Rw;
        S(i).Rt=sum(Ro)+sum(Rw);
    else
        %S(i).tr=[];
        %S(i).Ro=[];
        %S(i).Rg=[];
        %S(i).Rw=[];
        %S(i).Rt=0;
    end
    
end

% Plot for QC'ing.
figure(1); clf;
plot( -68.06,-38.96,'rd'); hold on;
for i=1:length(S)
    plot(S(i).Slon,S(i).Slat,'ob');
    plot(S(i).Wlon,S(i).Wlat,'-b');
end

% Save the data.
save('/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat','S');



