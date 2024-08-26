function [Vc]=getVdata(OP,Tv,type_flag)
  % Simple routine to get injected/extracted volumes from operational data structures.
  
  % Preallocate space for the output vector.
  Vc=zeros(size(Tv));
  
  % Loop over the operational data structure.
  for i=1:length(OP)
      
      % Get stimulation info for this one HF well.
      if(strcmpi(type_flag,'HF'))
          vt=interp1([OP(i).Ts;OP(i).Te+0.1;max(Tv)],[0;OP(i).Vt;OP(i).Vt],Tv,'linear',0);
          vt=[vt(1), diff(vt)];
          vt(vt<0)=0;
      elseif(strcmpi(type_flag,'INJ'))
          vt=interp1([0;OP(i).Tv],[0;cumsum(OP(i).Vw)],Tv,'linear',0);
          vt=[vt(1), diff(vt)];
          vt(vt<0)=0;
      elseif(strcmpi(type_flag,'PRD'))
          vt=interp1([0;OP(i).Tv],[0;cumsum(OP(i).Vw+OP(i).Vo+OP(i).Vg)],Tv,'linear',0);
          vt=[vt(1), diff(vt)];
          vt(vt<0)=0;
      end
      
      % Add this well into the whole dataset.
      Vc=Vc+vt;
  end
  
return