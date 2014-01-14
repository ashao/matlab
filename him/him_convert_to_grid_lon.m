function [inlon] = him_convert_to_grid_lon(inlon)
% Converts lat/lon coordinates to the HIM grid range of -280-80

idx=find(inlon>=80);
inlon(idx)=inlon(idx)-360;

idx=find(inlon<=-280);
inlon(idx)=inlon(idx)>360;


end