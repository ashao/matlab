function Dataset = nc_getvarinfo ( arg1, arg2 )
%NC_GETVARINFO  Returns metadata about a specific NetCDF variable.
%
%   VINFO = NC_GETVARINFO(NCFILE,VARNAME) returns a metadata structure 
%   about the variable VARNAME in the netCDF file NCFILE.
%
%   VINFO will have the following fields:
%
%       Name      - A string containing the name of the variable.
%       Datatype  - The datatype of the variable.
%       Unlimited - Either 1 if the variable has an unlimited dimension or 
%                   0 if not.
%       Dimension - a cell array with the names of the dimensions upon 
%                   which this variable depends.
%       Size      - Size of the variable.
%       Attribute - An array of structures corresponding to the attributes 
%                   defined for the specified variable.
%                         
%    Each "Attribute" element is a struct itself and contains the following 
%    fields.
%
%       Name      - A string containing the name of the attribute.
%       Datatype  - The datatype of the variable.
%       Value     - Value of the attribute.
%
%   See also nc_info.

if ~ischar(arg1)
    warning('SNCTOOLS:NC_GETVARINFO:deprecatedSyntax', ...
            'Using numeric IDs as arguments to NC_GETVARINFO is a deprecated syntax.');
end

backend = snc_read_backend(arg1);
switch(backend)
	case 'tmw'
		Dataset = nc_getvarinfo_tmw('',arg1,arg2);
	case 'java'
		Dataset = nc_getvarinfo_java(arg1,arg2);
	case 'mexnc'
		Dataset = nc_getvarinfo_mexnc(arg1,arg2);
    case 'tmw_hdf4'
        Dataset = nc_getvarinfo_hdf4(arg1,arg2);
    case 'tmw_hdf4_2011a'
        Dataset = nc_getvarinfo_hdf4_2011a(arg1,arg2);
	otherwise
		error('SNCTOOLS:nc_info:unhandledBackend', ...
		      '%s is not a recognized backend.', backend);
end



