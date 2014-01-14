function [ array ] = ncextractfield( ncfile,array,field,extent )
    %NCEXTRACTFIELD Extract the name of a field and given extent from a
    % Netcdf file
    %   Input:
    %       Required
    %       --------
    %       ncfile (string): path/filename of Netcdf file
    %       array (string): Structure to store data. If empty, new
    %                       structure is made
    %       field (string): name of field (e.g. 'mn_h' for offtrac layer
    %                       thickness)
    %       Optional
    %       --------
    %       extent (string): Data indices to extract (e.g. for offtrac,
    %                        (1,:) extracts all information in day 1 of
    %                        model output)
    
    if isempty(array)
        array=struct;
    end
    
    if nargin<4
        ncload(ncfile,field);
    else
        ncload2(ncfile,[field extent]);
    end
    
    array=setfield(array,field,eval(field));
    
    
end
