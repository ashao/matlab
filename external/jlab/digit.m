function[xdigit]=digit(x,n,flag,flag2)
%DIGIT  Returns the specified digit(s) of input numbers.
%
%   DIGIT(X,N) returns the Nth digit(s) of X. N=0 returns the ones digit, 
%   N=1 the tens digit, N=-1 the tenths digit, etc. X may be either a     
%   number or a vector, as may N. If X is a column (row) vector and       
%   length(N)>1, the N(1)th, N(2)th, ... digits are put in the 1st, 2nd,  
%   ... columns (rows) of the output.  
%          Example: DIGIT([123;456],[2 1])=[12;45];           
%                                                                         
%   DIGIT(X,N,FLAG) determines whether the output are strings             
%   (FLAG='str', the default) or numbers (FLAG='num').                    
%                                                                         
%   For string output, DIGIT(..., FLAG2) where FLAG2='spaces' gives empty 
%   digits filled with spaces (' 1'), while FLAG2='zeros' (the default)   
%   gives empty digits filled with zeros ('001'). 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
nin=n;
xin=x; 

btrans=0;
if size(x,1)==1 && size(x,2)>1
	x=x';
	btrans=1;
end

bstring=1;
bspace=0;
if nargin==3
	if strcmp(flag(1:3),'str')
		bstring=1;
	elseif strcmp(flag(1:3),'num')
		bstring=0;
	elseif strcmp(flag(1:3),'spa')
		bspace=1;
	elseif strcmp(flag(1:3),'zer')
		bspace=0;
	end
end

if nargin==4
	if strcmp(flag2(1:3),'str')
		bstring=1;
	elseif strcmp(flag2(1:3),'num')
		bstring=0;
	elseif strcmp(flag2(1:3),'spa')
		bspace=1;
	elseif strcmp(flag2(1:3),'zer')
		bspace=0;
	end
end

for i=1:length(nin)
	n=nin(i);
	x=row2col(xin);

	if n<1
		a=-n;
		x=x*10^a;
		n=0;
	end

	x=floor(x./(10^(n)));
	if bstring
		xdigit(:,i)=(int2str(x-floor(x./(10))*10))';
		if bspace && n~=0
			iii=find(xin<10^n)';
			if ~isempty(iii)
				xdigit(iii,i)=char(32+0*iii);
			end
		end
	else
		xdigit(:,i)=(x-floor(x./(10))*10);
	end
end

if btrans
	xdigit=xdigit';
end






function[col]=row2col(row)
%ROW2COL  Given any vector, returns it as a column vector.
 
if size(row,2)>size(row,1)
        if ~ischar(row)
                col=conj(row');
        else
                col=row';
        end
else
        col=row;
end

