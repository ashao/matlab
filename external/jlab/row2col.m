function[col]=row2col(row)
%ROW2COL  Given any vector, returns it as a column vector.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1999, 2004 J.M. Lilly --- type 'help jlab_license' for details      

if size(row,2)>size(row,1)
	if ~ischar(row)
     		col=conj(row');
	else 	
		col=row';
	end
else
	col=row;
end





