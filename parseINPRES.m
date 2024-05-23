function [lat,lon,dep,T,M]=parseINPRES(filename,lat_L,lon_L,dep_L,T_L,M_L)
  % Simple routine to parse csv files from INPRES.
  
  % Load in the EQ data.
  % yyyy, mm, dd, hh, mm, ss, lat, lon, dep, mag
  command=['cat ', filename, ' | ', ... % Load in all the EQ data files.
           'awk -F, ''/Source/{next} {print $2, $3, $4, $5, $6, $7 }'' | ' ... % Grab just the necessary fields.
           'awk     ''{split($1,d,""); split($2,t,":"); print d[1]d[2]d[3]d[4], d[5]d[6], d[7]d[8], t[1],t[2],t[3], $3, $4, $5, $6}'' > temp.EventParse']; % Split the time/date field into its respective parts.
  system(command);
  data=load('temp.EventParse');
  system('rm -f temp.EventParse');
  T=datetime(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6));
  lat=data(:,7);
  lon=data(:,8);
  dep=data(:,9);
  M=data(:,10);
  
  % Filter spatially (lateral).
  if(lat_L~=-1)
      I=inpolygon(lon,lat,lon_L,lat_L);
      T=T(I);
      M=M(I);
      lat=lat(I);
      lon=lon(I);
      dep=dep(I);
  end
  
  % Filter spatially (depth).
  if(dep_L~=-1)
      I=(dep>=min(dep_L))&(dep<=max(dep_L));
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
  
  % Sort chronologically.
  [~,I]=sort(T);
  T=T(I);
  M=M(I);
  lat=lat(I);
  lon=lon(I);
  dep=dep(I);
  
return


