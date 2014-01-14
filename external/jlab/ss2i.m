function[i]=ss2i(s1,s2)
%SS2I Convert from (s1,s2) notation to i=1:4 notation.
%
%   I=SS2I(S1,S2) where S1 and S2 are each +/- 1 returns a number I 
%   from 1 to 4 according to  
%
%        1 <==> (+,+)
%        2 <==> (-,-)
%        3 <==> (+,-)
%        4 <==> (-,+)
%  
%   See also I2SS.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    
  
 
  if s1==1  &&  s2==1
    i=1;
  elseif s1==-1  &&  s2==-1
    i=2;
  elseif s1==1  &&  s2==-1
    i=3;
  elseif s1==-1  &&  s2==1
    i=4;
  end
