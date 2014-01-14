function[str]=mat2table(x,prec)
%MAT2TABLE Converts a matrix of numbers into a LaTeX-style table.
%
%   STR=MAT2TABLE(MAT,PREC) converts a matrix of numbers into a string
%   matrix.  Columns of MAT are separated with ' & ' and the last
%   column in STR is followed by ' \\'.  This allows matrices to be
%   converted into a form that can be easily pasted into a LaTeX
%   document.
%
%   PREC specifies the precision of the last retained digit. For
%   example, PREC = -2 specifies that the hundredths digit is to be
%   retained, and so on.  PREC may either be a scalar, or an array
%   indicating the precision to use for each column of MAT.
%
%   See also VINDEX, DIGIT
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1999, 2004 J.M. Lilly --- type 'help jlab_license' for details        

if nargin==1
   prec=0;
end
if length(prec)==1
	prec=prec*ones(size(x(1,:)));
end
%put a \\ at end of all rows
str=char(real('\')*ones(size(x(:,1))));
str=[str str];


%work backwards
andv=char(real('&')*ones(size(x(:,1))));
blankv=char(real(' ')*ones(size(x(:,1))));
for i=size(x,2):-1:1
	if i==size(x,2)
	     	str=[vnum2str(x(:,i),prec(i),'spaces') blankv str];
	else
		str=[vnum2str(x(:,i),prec(i),'spaces') blankv andv blankv str];
	end
end



