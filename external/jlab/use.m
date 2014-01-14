function[]=use(x)
%USE  Copies structure fields into named variables in workspace.
%
%   USE STRUCT, where STRUCT is a structure, copies all fields of the
%   form STRUCT.X into variables named X.
%  
%   This is useful for handling multiple datasets with the same
%   variable names.  The structures can be then kept in memory and
%   'mapped' into variables as needed.
%
%   See also MAKE, MATSAVE.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2006 J.M. Lilly --- type 'help jlab_license' for details    

%'if ~exist (' x ')==1, load 'x'  

clear str
str{1}    =['if ~isempty(' x '),'];
str{end+1}=['  ZZFNAMES=fieldnames(' x ');' ];
str{end+1}='  for ZZi=1:length(ZZFNAMES),';
str{end+1}=[' 	  eval([ZZFNAMES{ZZi}, ''=getfield(' x ',ZZFNAMES{ZZi});'']);'];
str{end+1}='  end;';
str{end+1}='else;';
str{end+1}='  disp([''Contains no data.'']);'; 
str{end+1}='end;';
str{end+1}='clear ZZi ZZFNAMES';

str=strs2sray(str);
evalin('caller',str)



