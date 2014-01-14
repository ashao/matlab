function[bool,var]=isassignment(str)
%ISASSIGNMENT Checks if a string is a variable assignment, e.g. 'x=cos(t);'
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  

%isassignment(['x==10'])
%[n,var]=isassignment(['x=10'])
%isassignment(['x=10;                       ';'title(''Tests whether x=10''])'])
      
for i=1:size(str,1)
  index=findunquoted('=',str(i,:));
  
  %check that this is not part of a conditional expression
  bright=real(str(i,index+1))~=real('=');
  bleft= real(str(i,index-1))~=real('=') && ...
         real(str(i,index-1))~=real('~') && ...
         real(str(i,index-1))~=real('<') && ...
         real(str(i,index-1))~=real('>');
  index2=find(bleft.*bright);
  if ~isempty(index2)
    index=index(index2);
  else 
    index=[];
  end
  
  bool(i,1)=length(index);
  if bool(i,1)>0
    var{i}=str(i,1:index(1)-1);
  else
    var{i}=[];
  end
end


if length(var)==1
  var=var{1};
end
