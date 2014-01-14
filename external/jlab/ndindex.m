function[y]=ndindex(index,x,dim)
%NDINDEX  Indexes a multidimensional array along a specified dimension.
%
%   Y=NDINDEX(INDEX,X,DIM) indexes the multidimensional array X
%   along dimension DIM.  This is equivalent to
%		
%		    1 2       DIM     DIMS(X)
%		    | |        |         |
%		Y=X(:,:, ... INDEX, ..., :);		
%
%   where the location of INDEX is specified by DIM.
%
%   See also SQUEEZE, DIMS, PERMUTE, SHIFTDIM.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details        
  

%You would think Matlab would provide a simpler way to do this.

str='y=x(';
ndx=length(find(size(x)>1));
for i=1:ndx
    if i~=dim
        str=[str ':,'];
    else
	str=[str 'index,'];
    end
end
str=[str(1:end-1) ');'];
eval(str);
	
