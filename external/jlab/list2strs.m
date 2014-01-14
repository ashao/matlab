function[y]=list2strs(x)
%LIST2STRS   Converts a comma-delimited list into a cell array of strings.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        
y=[];

index=find(real(x)==real(','));

if ~isempty(index)
  a=[1 index+1];
  b=[index-1 length(x)];
  for i=1:length(a)
    y{i}=deblank(x(a(i):b(i)));
  end
end

y=packstrs(y);




if 0
    if ~isempty(index)
        if index(end)==length(x)
            index=index(1:end-1);
            x=x(1:end-1);
        end
    end
    if ~isempty(index)
        if index(1)==1
            index=index(2:end);
        end
    end
end
