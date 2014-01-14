function[x]=flushright(x)
%FLUSHRIGHT   Makes a blank-padded string matrix flush on the right
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details  
  
if strcmp(x,'--t')
    return
end

x=fliplr(flushleft(fliplr(x)));
