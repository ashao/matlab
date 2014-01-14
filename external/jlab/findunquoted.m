function[index]=findunquoted(s1,s2)
%FINDUNQUOTED Finds unquoted instances of one string inside another
%
%   FINDUNQUOTED(A,B) where A and B are strings, finds unquoted
%   instances of A inside B, and returns and an array specifiying 
%   the position within B of the first letter of each occurrence.
%
%   For example
%     findunquoted('apple',['banana apple ''apple'' apple'])
%   returns [8 22]. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
%findunquoted('==',['plot(1:100);x==10;title([''Tests whether x==10''])'])  

index=strfind(s1,s2);
bool=isquoted(s1,s2);
index=index(~bool);

