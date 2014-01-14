function[varargout]=offset(varargin)
%OFFSET Offsets matrix rows or columns
%
%   MATOUT = OFFSET(MATIN,SHIFT,DIM), where SHIFT is a scalar, offsets
%   dimension DIM of MATIN by SHIFT. 
%  
%   Use OFFSET to produce 'waterfall' plots of data. ex:                  
%   [x,y,z]=peaks; plot(offset(z,2,2))                          
%
%   OFFSET(MATIN,SHIFT,DIM) where MATIN is a named matrix
%   overwrites it with MATIN.
%
%   See also YOFFSET, XOFFSET  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
%  JML Make this work for N-D

dim=varargin{end};
offset=varargin{end-1};

for i=1:length(varargin)-2
   varargout{i}=offset1(varargin{i},offset,dim);
end

eval(to_overwrite(nargin-2))
    
function[y]=offset1(x,offset,dim)

N=size(x,dim);

if isscalar(offset)
  offset=offset*(0:N-1);
end

if N~=length(offset)
  error('Length of offset vector must equal SIZE(X,DIM)')
end

%Reverse dim for calling VADD
if dim==2
  dim=1;
elseif dim==1
  dim=2;
end


%make sure they are oriented the right way
if dim==1
  offset=conj(offset(:)');
elseif dim==2
  offset=offset(:);
end
y=vadd(x,offset,dim);

       
       
       
