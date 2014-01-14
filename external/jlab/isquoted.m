function[bool,n]=isquoted(str1,str2)
%ISQUOTED Checks to see if a string is quoted in a longer expression
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details
  
%isquoted('==',['x==10'])
%isquoted('==',['title([''Tests whether x==10''])'])

str2=str2(:)';
index=strfind(str1,str2);
n=length(index);
if isempty(index)
  bool=0;
else
  quotebool=real(str2)==real(''''); % 1 for occurences of a quoteation
  
  for i=1:length(index)
    leftquotes=sum(quotebool(1:index(i)));
    rightquotes=sum(quotebool(index(i):end));

    %if something is quoted, we'll see an odd number of quotes each way

    if isodd(leftquotes) && isodd(rightquotes)
      bool(i,1)=1;
    else 
      bool(i,1)=0;
    end
  end
end
