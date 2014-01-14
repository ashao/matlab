function values = snc_pp_strings(jobj,jdata,shape)
% Post process NC_STRING data into cell arrays.
if isempty(jdata)
    values = {''};
    return;
elseif strcmp(version('-release'),'14')
    % In R14, we must use the 2.2.x release of java.  No access to the
    % "getObject" method.  Assuming a single-valued string.
    values = {jobj.getStringValue()};
    return;
end


% Java says that the variable is laid out in row-major order.
if numel(shape) == 1
    values = cell([1 shape]);
else
    values = cell(shape);
end

for j = 1:prod(shape)
    values{j} = jdata.getObject(j-1);
end
