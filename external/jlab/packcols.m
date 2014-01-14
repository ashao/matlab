function[h]=packcols(i1,i2,i3)
%PACKCOLS   Squeeze together all subplot columns of the current figure.
%  
%   PACKCOLS(M,N) squeezes together all columns in the current figure,
%   which has M rows and N columns generated with SUBPLOT. This is
%   used when all subplots in a given column share a common y-axis.
%   Y-axis tick labels for the interior subplots are removed.
%
%   H=PACKCOLS(M,N) returns a vector of handles H into the subplots.   
%
%   PACKCOLS(H,M,N) also works, where H is a vector of handles.
%
%   After calling PACKCOLS, do not call SUBPLOT again, as this will
%   destroy the subplots; instead, access the subplots through
%   AXES(H(I)) where I is an integer. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2008 J.M. Lilly --- type 'help jlab_license' for details        

%Depricated
%    If the
%   XTICKMODE property is set to manual (i.e. if the XTICKS have been
%   explicitly specified) then the last value of XTICKS is set to
%   blanks in all columns save the rightmost.
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

if n>1
  
  for i=2:n
    for j=1:m
      index=(j-1)*n+i;
      axes(h(index))
      noylabels,ylabel('')
    end
  end

  axes(h(1))
  pos1=get(gca,'position');
     
  axes(h(end))
  pos2=get(gca,'position');
      
  dx=pos2(1)+pos2(3)-pos1(1);   %total width to work with 
  
  xwidth=frac(dx-dd*(n-1),n);  %This is how wide each axis should be 
  xleft=zeros(n,1);
  xleft(1)=pos1(1);
  
  for i=2:n;
      xleft(i)=xleft(i-1)+xwidth+dd;
  end
  %xleft=flipud(xleft);
  
  
  for i=1:n
    for j=1:m
      index=(j-1)*n+i;
     
      axes(h(index))
      pos=get(gca,'position');
      isboxoff=strcmp(get(gca,'box'),'off');
      
      pos(1)=xleft(i);
      pos(3)=xwidth;
     
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
