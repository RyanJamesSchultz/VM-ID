function [S]=SAF(OP,latEQ,lonEQ,tEQ,dR,Pr,dT,Pt,decay_flag,assoc_flag)
  % A function that performs a rudimentary spatiotemporal association.
  
  % Predefine the output structure.
  S=struct('Ieq',[],'Ihf',[],'IDhf',[],'Ir',[],'It',[],'L',[],'Lr',[],'Lt',[],'Req',[],'Twe',[],'Ns',[],'Nt',[],'N',[],'S',[]);
  
  % Loop over all of the operations and find associated earthquakes.
  for j=1:length(OP)
      
      % Pull out the requisite information.
      latOP=OP(j).Slat;
      lonOP=OP(j).Slon;
      tsOP=OP(j).T(1);
      teOP=OP(j).T(end);
      
      % Compute epicentral distances from surface location to EQ (km).
      Req=Geoid_Distance(latOP,lonOP, latEQ,lonEQ, 'elliptical')*6371*pi()/180; % km.
      Teq=daysdif(teOP,tEQ); % days from shut-in.
      %dTl=daysdif(tsOP,teOP)+dT;
      
      % Filter spatiotemporally.
      if(strcmpi(decay_flag,'exp'))
          Jr=Req<=dR(1)+7*dR(2);
          Jt=(tEQ>=tsOP)&(Teq<=(dT(1)+6*dT(2)));
          J=Jt&Jr;
      elseif(strcmpi(decay_flag,'exp2'))
          Jr=Req<=dR(1)+5*dR(2);
          Jt=(tEQ>=tsOP)&(Teq<=(dT(1)+6*dT(2)));
          J=Jt&Jr;
      elseif(strcmpi(decay_flag,'power'))
          Jr=Req<=Inf;
          Jt=(tEQ>=tsOP)&(Teq<=Inf);
          J=Jt&Jr;
      else
          Jr=Req<=dR(2);
          Jt=(tEQ>=tsOP)&(Teq<=dT(2));
          J=Jt&Jr;
      end
      
      % Compute SAF-scores for each event.
      Lr=Lscore(Req(J),Pr,dR(1),dR(2),decay_flag);
      Lt=Lscore(Teq(J),Pt,dT(1),dT(2),decay_flag);
      L=Lr.*Lt;
      
      % Stuff results into the output structure.
      S(j).Ihf=j;
      S(j).IDhf=OP(j).IDe;
      S(j).Ieq=find(J);
      S(j).Ir=find(Jr);
      S(j).It=find(Jt);
      S(j).L=L;
      S(j).Lr=Lr;
      S(j).Lt=Lt;
      S(j).Req=Req(J);
      S(j).Twe=Teq(J);
      S(j).N=sum(J);
      S(j).Ns=sum(Jr);
      S(j).Nt=sum(Jt);
      S(j).S=sum(L);
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
          J=[]; Lt=[];
          for j=1:length(S)
              I=S(j).Ieq==i;
              if(any(I))
                  J=[J,j]; Lt=[Lt,S(j).L(I)];
              end
          end
          
          % Find the well that is best associated.
          [~,jb]=max(Lt); Jb=J(jb); J(jb)=[];
          
          % Remove the event from the all of the 'runner up' wells.
          for j=J
              I=S(j).Ieq==i;
              
              % Remove the event
              S(j).Ieq(I)=[];
              S(j).Ir(S(j).Ir==i)=[];
              S(j).It(S(j).It==i)=[];
              S(j).L(I)=[];
              S(j).Lr(I)=[];
              S(j).Lt(I)=[];
              S(j).Req(I)=[];
              S(j).Twe(I)=[];
              
              % Update some parameters.
              S(j).N =length(S(j).Ieq);
              S(j).Ns=length(S(j).Ir);
              S(j).Nt=length(S(j).It);
              S(j).S =sum(S(j).L);
          end
          
      end
      
  end
  
return




% Individual event-well spatial/temporal-score sub-routine.
function [P]=Lscore(x,p1,x1,x2,decay_flag)
  % A function to calculate a spatial-score or temporal-score.
  
  % Check the user-input for functional form of score decay.
  if(strcmpi(decay_flag,'exp'))
      
      % Compute the score as piecewise: linear over first interval, exponential taper over second.
      p1i=1-p1;
      P=1-(p1i/x1)*x; % scores of 1-p1 between 0-x1.
      P(x>x1)=p1*exp(-(x(x>x1)-x1)/x2); % scores of p1-0 between x1-x2.
      P(x<=0)=1;
      
  elseif(strcmpi(decay_flag,'exp2'))
      
      % Compute the score as piecewise: linear over first interval, Gaussian taper over second.
      p1i=1-p1;
      P=1-(p1i/x1)*x; % scores of 1-p1 between 0-x1.
      P(x>x1)=p1*exp(-(x(x>x1)-x1).^2/(2*x2^2)); % scores of p1-0 between x1-x2.
      P(x<=0)=1;
      
  elseif(strcmpi(decay_flag,'power'))
      
      % Compute the score as piecewise: linear over first interval, power taper over second.
      P(x<=x1)=1; % scores of 1 between 0-x1.
      P(x>x1)=p1*x(x>x1).^x2; % scores of p1-0 between x1-x2.
      
  else
      
      % Compute the score as a bilinear function.
      p1i=1-p1;
      P=1-(p1i/x1)*x; % scores of 1-p1 between 0-x1.
      P(x>x1)=p1-(p1/(x2-x1))*(x(x>x1)-x1); % scores of p1-0 between x1-x2.
      P(x<=0)=1;
      P(x>x2)=0;
  end
  
return


