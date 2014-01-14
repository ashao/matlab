function attribute = nc_getattsinfo_tmw(ncfile,ncid,varid,attnum)
% NC_GET_ATTRIBUTE_STRUCT_TMW:  Returns a NetCDF attribute as a structure
%
% You don't want to be calling this routine directly.  Just don't use 
% it.  Use nc_attget instead.  Go away.  Nothing to see here, folks.  
% Move along, move along.


attribute = struct('Name','','Nctype','','Datatype','','Value',NaN);


attname = netcdf.inqAttName(ncid, varid, attnum);
attribute.Name = attname;

[att_datatype] = netcdf.inqAtt(ncid, varid, attname);
attribute.Nctype = att_datatype;
switch(att_datatype)
    case nc_nat
        attribute.Datatype = '';
    case nc_byte
        attribute.Datatype = 'int8';
    case nc_ubyte
        attribute.Datatype = 'uint8';
    case nc_char
        attribute.Datatype = 'char';
    case nc_short
        attribute.Datatype = 'int16';
    case nc_ushort
        attribute.Datatype = 'uint16';
    case nc_int
        attribute.Datatype = 'int32';
    case nc_uint
        attribute.Datatype = 'uint32';
    case nc_int64
        attribute.Datatype = 'int64';
    case nc_uint64
        attribute.Datatype = 'uint64';
    case nc_float
        attribute.Datatype = 'single';
    case nc_double
        attribute.Datatype = 'double';
    case 12
        attribute.Datatype = 'string';
    otherwise
        attribute.Datatype = '';
end

switch att_datatype
    case 0
        attribute.Value = NaN;
    case { nc_char, nc_int64, nc_uint64 }
        attribute.Value=netcdf.getAtt(ncid,varid,attname);
    case { nc_double, nc_float, nc_int, nc_short, nc_byte, nc_ubyte, nc_ushort, nc_uint }
        attribute.Value=netcdf.getAtt(ncid,varid,attname,'double');
    case 12 % string
        attribute = nc_getattsinfo_java(ncfile,ncid,varid,attname);

    otherwise
        warning ( 'SNCTOOLS:nc_getattsinfo:tmw:unhandledDatatype', ...
            'The datatype for attribute ''%s'' (%d) is not currently handled by SNCTOOLS.', ...
            attname, att_datatype );
        attribute.Value = [];
end


return


%---------------------------------------------------------------------------
function info = nc_getattsinfo_java(ncfile,ncid,varid,attname)

import ucar.nc2.dods.*     
import ucar.nc2.*  

jncid = NetcdfFile.open(ncfile);

group_name = netcdf.inqGrpNameFull(ncid);

if strcmp(group_name,'/')
    gid = jncid;
else
    group_name = group_name(2:end);
    
    gid = jncid;
    
    % Strip off the leading slash to make REGEXP's life a bit easier.
    groups = regexp(group_name,'/','split');
    for j = 1:numel(groups)
        gid = gid.findGroup(groups{j});
    end
    
end


if isa(gid,'ucar.nc2.NetcdfFile') && (varid == -1)
    jatt = gid.findGlobalAttribute(attname);
elseif (varid == -1)
    jatt = gid.findAttribute(attname);   
else
    varname = netcdf.inqVar(ncid,varid);
    
    jvarid = gid.findVariable(varname);
    jatt = jvarid.findAttribute(attname);
end

info = nc_getattinfo_java(jatt);
close(jncid);

