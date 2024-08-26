% Script that will parse VM-HF well header files and injection files into
% my Matlab structure.
clear;

% Get the well header information.
Th=readtable('/Users/rschultz/Desktop/VacaMuerta/data/newdata/listado-de-pozos-cargados-por-empresas-operadoras.csv');

% Get the directional-survey and stimulation data.
Td=[]; Ti=[]; Tp=[];
%Td=readtable('DirSurv.csv');
Ti=readtable('/Users/rschultz/Desktop/VacaMuerta/data/newdata/datos-de-fractura-de-pozos-de-hidrocarburos-adjunto-iv-actualizacin-diaria.csv');
%Tp=readtable('Prod2.csv');


% Define the data strucutre.
S=struct('UWI',[],'IDe',[],'type',[],'Slat',[],'Slon',[],'Lh',[],'Ns',[],'T',[],'V',[],'P',[],'Vt',[],'Wlat',[],'Wlon',[],'Wtvd',[]);

% Loop over all of the well IDs.
for i =1:size(Th,1)
    
    % Stuff the header information into the injection data structure.
    uwi=Th{i,2};
    type=Th{i,47};
    id=Th{i,1};
    S(i).IDe=id;
    S(i).UWI=uwi{1};
    S(i).type=type{1};
    S(i).Slat=Th{i,12};
    S(i).Slon=Th{i,11};
    
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
    Ii=find(strcmpi(Ti{:,3},S(i).UWI));
    if(~isempty(find(Ii,1)))
        T=[Ti{Ii(1),19} Ti{Ii(1),20}];
        P=Ti{Ii(1),17}*0.0068947573; % psi to MPa
        V=[0 Ti{Ii(1),15}]; % m³.
        V(isnan(V))=0;
        S(i).T=T;
        S(i).V=V;
        S(i).P=P;
        S(i).Vt=Ti{Ii(1),15}; % m³.
        S(i).Lh=Ti{Ii(1),10};
        S(i).Ns=Ti{Ii(1),11};
    else
        S(i).T=[];
        S(i).V=[];  % m³.
        S(i).P=[];  % MPa
        S(i).Vt=[]; % m³.
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

% Remove the non-HF (or missing) data.
i=1;
while(i<=length(S))
    if(isempty(S(i).T))
        S(i)=[];
    elseif((S(i).Slat<-40)||(S(i).Slat>-37))
        S(i)=[];
    elseif(isempty(S(i).Vt))
        S(i)=[];
    else
        if(S(i).Vt==0)
            S(i)=[];
        else
            i=i+1;
        end
    end
end

% Deal with non-uniqueness in the header database.
[~,I]=unique({S.UWI});
S=S(I);

% Plot for QC'ing.
figure(1); clf;
plot( -68.06,-38.96,'rd'); hold on;
for i=1:length(S)
    plot(S(i).Slon,S(i).Slat,'ob');
    plot(S(i).Wlon,S(i).Wlat,'-b');
end

% Save the data.
save('/Users/rschultz/Desktop/VacaMuerta/data/INJ2.mat','S');



