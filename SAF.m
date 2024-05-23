function [S]=SAF(OP,latEQ,lonEQ,tEQ,dR,dT,assoc_flag)
  % A function that performs a rudimentary spatiotemporal association.
  
  % Predefine the output structure.
  S=struct('Ieq',[],'Ihf',[],'Is',[],'It',[],'P',[],'Ps',[],'Pt',[],'Req',[],'Twe',[],'Ns',[],'Nt',[],'N',[],'S',[]);
  
  % Loop over all of the operations and find associated earthquakes.
  for i=1:length(OP)
      
      % Pull out the requisite information.
      latOP=OP(i).Slat;
      lonOP=OP(i).Slon;
      tsOP=OP(i).T(1);
      teOP=OP(i).T(end);
      
      % Compute epicentral distances from surface location to EQ (km).
      Req=Geoid_Distance(latOP,lonOP, latEQ,lonEQ, 'elliptical')*6371*pi()/180; % km.
      Teq=daysdif(tsOP,tEQ); % days.
      dTl=daysdif(tsOP,teOP)+dT;
      
      % Filter spatiotemporally.
      Ir=Req<=dR(2);
      It=(Teq>=0)&(Teq<=dTl);
      I=It&Ir;
      
      % Compute SAF-scores for each event.
      Pr=Pscore(Req(I),0.95,dR(1),dR(2));
      Pt=Pscore(Teq(I),1.00,dTl-dT,dTl);
      P=Pr.*Pt;
      
      % Stuff results into the output structure.
      S(i).Ihf=OP(i).IDe;
      S(i).Ieq=find(I);
      S(i).Is=find(Ir);
      S(i).It=find(It);
      S(i).P=P;
      S(i).Ps=Pr;
      S(i).Pt=Pt;
      S(i).Req=Req(I);
      S(i).Twe=tEQ(I)-teOP;
      S(i).N=sum(I);
      S(i).Ns=sum(Ir);
      S(i).Nt=sum(It);
      S(i).S=sum(P);
  end
  
  % Keep only the unique associations, if flagged to.
  if(strcmpi(assoc_flag,'unique'))
      
      % Loop over each earthquake in the catalogue.
      for i=1:length(latEQ)
          
          % Check if the event has multiple associations, ignore if it's already unique or unassociated.
          I=vertcat(S.Ieq);
          if(length(find(I==i))<=1)
              continue;
          end
          
          % Find all of the wells that are associated.
          J=[]; Pt=[];
          for j=1:length(S)
              I=S(j).Ieq==i;
              if(any(I))
                  J=[J,j]; Pt=[Pt,S(j).P(I)];
              end
          end
          
          % Find the well that is best associated.
          [~,jb]=max(Pt); Jb=J(jb); J(jb)=[];
          
          % Remove the event from the all of the 'runner up' wells.
          for j=J
              I=S(j).Ieq==i;
              
              % Remove the event
              S(j).Ieq(I)=[];
              S(j).Is(S(j).Is==i)=[];
              S(j).It(S(j).It==i)=[];
              S(j).P(I)=[];
              S(j).Ps(I)=[];
              S(j).Pt(I)=[];
              S(j).Req(I)=[];
              S(j).Twe(I)=[];
              
              % Update some parameters.
              S(j).N =length(S(j).Ieq);
              S(j).Ns=length(S(j).Is);
              S(j).Nt=length(S(j).It);
              S(j).S =sum(S(j).P);
          end
          
      end
      
  end
  
return




% SAF-score sub-routine.
function [P]=Pscore(x,p1,x1,x2)
  % A function to calculate a spatiotemporal association score.
  
  % Compute the score as a bilinear function.
  p1i=1-p1;
  P=1-(p1i/x1)*x; % scores of 1-p1 between 0-x1.
  P(x>x1)=p1-(p1/(x2-x1))*(x(x>x1)-x1); % scores of p1-0 between x1-x2.
  P(x<0)=0;
  P(x>x2)=0;
  
return


