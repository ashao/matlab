function [ idx ] = sel_col(prop,name)

    [m n] = size(prop);
    for ii=1:m
       prop(ii,:);
        if strcmpi(strtrim(prop(ii,:)),strtrim(name))
            idx=ii;
            break
        end
%         error('Cannot find %s\n',name);
        idx=[];
    end
    
    
end