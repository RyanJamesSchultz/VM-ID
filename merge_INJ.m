% Script that merges two (overlapping) VM-HF databases into one.
clear;

% Load in the two databases.
load('/Users/rschultz/Desktop/VacaMuerta/data/INJ1.mat','S'); S1=S;
load('/Users/rschultz/Desktop/VacaMuerta/data/INJ2.mat','S');

% Loop over all of the HF wells in the S1 database.
IDe=max([S.IDe])+100;
for i=1:length(S1)

    % Check if this HF well is in the S2 database.
    I=strcmpi({S.UWI},S1(i).UWI);
    if(any(I))    % This well is already in the S2 database.
       
        % Overwrite conflicting details.
        %S(I).T(1)=min([S1(i).T,S(I).T]);
        %S(I).T(2)=max([S1(i).T,S(I).T]);
        %S(I).Vt=max([S1(i).Vt,S(I).Vt]);
        %S(I).V(2)=S(I).Vt;
        %S(I).P=max([S1(i).P,S(I).P]);
        S(I).T=S1(i).T;
        S(I).V=S1(i).V;
        S(I).P=S1(i).P;
        S(I).Vt=S1(i).Vt;
        
        % Optionall report details.
        %S1(i).UWI
        %{S(I).UWI}
        %[S1(i).T;S(I).T]
        %[S1(i).Slat;S(I).Slat]
        %[S1(i).Slon;S(I).Slon]
        
    else          % This well is not in the S2 database yet.
        
        % Add it to the database and then iterate.
        S(end+1)=S1(i);
        S(end).IDe=IDe;
        IDe=IDe+1;
        
        % Optionally report details.
        %disp('Add')
        %S1(i).UWI
    end
end

% Save the data.
save('/Users/rschultz/Desktop/VacaMuerta/data/INJ.mat','S');

