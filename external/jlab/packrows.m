function[h]=packrows(i1,i2,i3)
%PACKROWS   Squeeze together all subplot rows of the current figure.
%  
%   PACKROWS(M,N) squeezes together all rows in the current figure,
%   which has M rows and N columns generated with SUBPLOT. This is
%   used when all subplots in a given column share a common x-axis.
%   X-axis tick labels for the interior subplots are removed.  
%
%   H=PACKROWS(M,N) returns a vector of handles H into the subplots.   
%
%   PACKROWS(H,M,N) also works, where H is a vector of handles.
%
%   After calling PACKROWS, do not call SUBPLOT again, as this will
%   destroy the subplots; instead, access the subplots through
%   AXES(H(I)) where I is an integer. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2008 J.M. Lilly --- type 'help jlab_license' for details        

%Depricated
%If the
%   YTICKMODE property is set to manual (i.e. if YTICK has been
%   explicitly specified) then the largest value of YTICK is set to
%   blanks in all rows save the uppermost.  

dd=.01;

if nargin==2
  m=i1;
  n=i2;
  h=subplots(m,n);
elseif nargin==3
  h=i1;
  m=i2;
  n=i3;
end


if m>1
  
  for i=1:n
    for j=1:m-1
      index=(j-1)*n+i;
      axes(h(index))
      noxlabels,xlabel('')
    end
  end

  axes(h(1))
  pos1=get(gca,'position');
     
  axes(h(end))
  pos2=get(gca,'position');
      
  dy=pos1(2)+pos1(4)-pos2(2);    %total height to work with 
  
  yheight=frac(dy-dd*(m-1),m);  %This is how tall each axis should be 
  ybottom=zeros(m,1);
  ybottom(1)=pos2(2);
  
  for j=2:m;
      ybottom(j)=ybottom(j-1)+yheight+dd;
  end
  ybottom=flipud(ybottom);
  
  
  for i=1:n
    for j=1:m
      index=(j-1)*n+i;
     
      axes(h(index))
      pos=get(gca,'position');
      isboxoff=strcmp(get(gca,'box'),'off');
      
      pos(2)=ybottom(j);
      pos(4)=yheight;
     
      set(gca,'position',pos,'box','on')
      if isboxoff
          boxoff
      end
    end
  end
end

if nargout==0
  clear h
end

