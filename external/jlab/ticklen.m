function[h]=ticklen(n)
%TICKLEN  Sets tick length of current axis.
%
%   TICKLEN(N) sets tick length of current axis to N times the default.
%  
%   TICKLEN([A B]) sets tick length of current axis to [A B].
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details

if strcmp(n,'--t')
    return
end

if length(n)==1
  set(gca,'ticklen',n*[0.0100    0.0250])
else
  set(gca,'ticklen',n)
end
