function [ idx,val ] = findrange( dat, range  )
%	FINDRANGE Returns data and indices in range

minval=min(range);
maxval=max(range);
idx=find( dat>=minval & dat<=maxval );
val=dat(idx);
end
