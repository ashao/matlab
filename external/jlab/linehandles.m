function[h,boolline,boolpatch]=linehandles(axh)
%LINEHANDES  Finds all line and patch handles from a given set of axes.
%
%   H=LINEHANDLES returns handles to all lines and patches associated
%   with the current axes.
%	 
%   H=LINEHANDLES(AX) returns handles to all lines and patches associated 
%   with the set of axes whose handle is AX.
%
%   [H,BLINE,BPATCH]=LINEHANDLES also returns two logical arrays of the
%   same size as H.  BLINE is true if the corresponding element of H is 
%   a line object, while BPATCH is true for patch objects.
%
%   See also PATCHHANDLES, AXESHANDLES
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2009 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==0
   axh=gca;
end
h=get(axh,'children');

hchild=[];
for i=1:length(h)
    if strcmp(get(h(i),'type'),'hggroup')
        hchild{i,1}=get(h(i),'children');
    end
end

if ~isempty(hchild)
    h=[h;vcellcat(hchild)];
end


htype=get(h,'type');
if ~iscell(htype)
    htypetemp=htype;
    clear htype
    htype{1}=htypetemp;
end

boolline=false(size(h));
boolpatch=false(size(h));
for j=1:length(h)
    boolline(j)= strcmp(htype{j},'line');
    boolpatch(j)= strcmp(htype{j},'patch');   
end

vindex(h,boolline,boolpatch,boolline|boolpatch,1);

