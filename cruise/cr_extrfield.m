function [ cruise ] = cr_extrfield( cruise, field )
%EXTRACTCRUISE Extracts specified field from cruise.data as separate field
%cruise.<field>

fieldid=cr_findfield(cruise,field);
field=field(isstrprop(field,'alphanum'));
cruise=setfield(cruise,field,cruise.data(:,fieldid));
end
