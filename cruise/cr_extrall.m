function [ cruise ] = cr_extrall( cruise )
%CR_EXTRALL Extracts all fields in the input array
% Input: cruise (struct) with fields data and prop
% Output: cruise (struct) with all fields listed in prop

for i=1:length(cruise.prop)
field=strtrim(lower(cruise.prop(i,:)));
cruise=cr_extrfield(cruise,field);
end

