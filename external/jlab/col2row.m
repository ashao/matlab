function[row]=col2row(col)
%COL2ROW  Given any vector, returns it as a row vector.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details

if size(col,1)>size(col,2)
	if ~ischar(col)
     		row=conj(col');
	else 	
		row=col';
	end
else
	row=col;
end



