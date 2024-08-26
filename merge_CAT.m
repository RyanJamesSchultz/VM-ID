% Script that merges two (overlapping) catalogues into one.
clear;

% Threshold values.
Tx=10;
Tt=Tx/6;

% Load the Neuquén Basin boundaries.
data=load('/Users/rschultz/Desktop/VacaMuerta/data/shapes_vaca_muerta/Neuquen.txt');
latB=data(:,2); lonB=data(:,1);

% Boundary polygons and stuff.
depB=-1;
tB=-1;
mB=-1;

% File names.
catFILE1='/Users/rschultz/Desktop/VacaMuerta/data/20240314_VM_quakes_all_soruces.csv';
catFILE2='/Users/rschultz/Desktop/VacaMuerta/data/final_base_phd.csv';

% Get the EQ catalogue data.
[EQlat1,EQlon1,EQdep1,EQtime1,EQmag1]=parseINPRES(catFILE1,latB,lonB,depB,tB,mB);
[EQlat2,EQlon2,EQdep2,EQtime2,EQmag2]=parseCorrea(catFILE2,latB,lonB,depB,tB,mB);

% Merge the two catalogues into one.
EQlat= [EQlat1; EQlat2];
EQlon= [EQlon1; EQlon2];
EQdep= [EQdep1; EQdep2];
EQtime=[EQtime1;EQtime2];
EQmag= [EQmag1; EQmag2];
EQids=[ones(size(EQmag1));2*ones(size(EQmag2))];

% Sort chronologically.
[~,I]=sort(EQtime);
EQtime=EQtime(I);
EQmag=EQmag(I);
EQlat=EQlat(I);
EQlon=EQlon(I);
EQdep=EQdep(I);
EQids=EQids(I);

% Remove list.
Ir=[];

% Check for redundancies within the merged catalogue.
for i=1:length(EQlat)
    
    % Find all of the (later) events within 2 seconds of this event.
    dt=seconds(EQtime-EQtime(i));
    I=find((dt<=Tt)&(dt>0));
    
    % Examine these temporally close events.
    if(~isempty(I))
        
        % Check if these candidate events are also spatially close.
        dR=Geoid_Distance(EQlat(i),EQlon(i),EQlat(I),EQlon(I),'elliptical')*6371*pi()/180;
        I(dR>Tx)=[];

        % Examine these spatiotemporally close events.
        if(~isempty(I))

            % Optionally report values to the screen.
            %i
            %I
            %[EQlat(i);EQlat(I)]
            %[EQlon(i);EQlon(I)]
            %[EQdep(i);EQdep(I)]
            %[EQtime(i);EQtime(I)]
            %[EQmag(i);EQmag(I)]
            %[EQids(i);EQids(I)]
            
            % Preferrentially keep the event from Correa-Otto et al., 2024.
            if(EQids(i)==2)
                Ir=[Ir,I];
            elseif(EQids(I)==2)
                Ir=[Ir,i];
            else
                Ir=[Ir,I];
            end
            
        end
    end

end

% Remove redundant events.
EQtime(Ir)=[];
EQmag(Ir)=[];
EQlat(Ir)=[];
EQlon(Ir)=[];
EQdep(Ir)=[];
EQids(Ir)=[];

% Save the data.
save('/Users/rschultz/Desktop/VacaMuerta/data/CAT.mat','EQlat','EQlon','EQdep','EQtime','EQmag');

% Plot.
figure(1); clf;
subplot(311);
plot(EQtime,EQmag,'o');
xlabel('Time'); ylabel('Magnitude');
subplot(3,1,[2 3]);
plot(EQlon,EQlat,'o');
xlabel('Longitude'); ylabel('Latitude');
