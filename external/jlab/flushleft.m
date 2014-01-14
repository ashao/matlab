function[x]=flushleft(x)
%FLUSHLEFT   Makes a blank-padded string matrix flush on the left
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details  
  
% Find orientation

if strcmp(x,'--t')
    return
end

for i=1:size(x,1)
  index=min(find(real(x(i,:))~=real(' ')));
  if ~isempty(index)
    x(i,:)=vshift(x(i,:),index-1,2);
  end
end
