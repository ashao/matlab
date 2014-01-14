function[s1,s2]=i2ss(i)
%I2SS Convert from i=1:4 notation to (s1,s2) notation.
%  
%   [S1,S2]=I2SS(I) where I is a number from 1 to 4 returns S1 and S2
%   each as +/- 1 according to   
%
%        1 <==> (+,+)
%        2 <==> (-,-)
%        3 <==> (+,-)
%        4 <==> (-,+)
%  
%   See also SS2I.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    
  
   

  switch i
   case 1 
    s1=1;
    s2=1;
   case 2 
    s1=-1;
    s2=-1;
   case 3
    s1=1;
    s2=-1;
   case 4
    s1=-1;
    s2=1;
  end
  

  

