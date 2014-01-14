function [ idx ] = cr_findfield( array, field )
    %FINDFIELD Returns column index of the given field in the array, returns
    %    null if nonexistent

    idx=1;
    field=strtrim(field);
    while ~strcmpi(strtrim(array.prop(idx,:)),field) && idx<length(array.prop)
        idx=idx+1;
    end

    if ~strcmpi(strtrim(array.prop(idx,:)),field);
        idx=[];
    end

end
