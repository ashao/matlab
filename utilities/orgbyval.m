function [ outtab uniqueval ] = orgbyval( orgvals, data )
% ORGBYVAL Constructs a data table whose columns all share the same value in orgval
% Input:
%	orgvals: values to sortby (length m)
%	data: vector of data values (length m)
% Output:
%	outtab: array whose columns are organized by orgvals
%	uniqueval: value for each column

uniqueval=unique(orgvals);
for i=1:length(uniqueval)
idx=find(orgvals==uniqueval(i));
outtab(1:length(idx),i)=data(idx);
end

end
