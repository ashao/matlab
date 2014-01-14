function nc_attput_tmw ( ncfile, varname, attribute_name, attval )
%NC_ATTPUT_TMW private function for writing attribute with TMW backend.

if isnumeric(ncfile) && isnumeric(varname)
    % the variable and file are already open, write the attribute now.
    write_attribute(ncfile,varname,attribute_name,attval);
    return
end

ncid  =netcdf.open(ncfile, nc_write_mode );

try
    library_version = netcdf.inqLibVers();
    library_version = str2double(library_version(1));
    
    netcdf.reDef(ncid);

    % If netcdf-4, make sure that the user did not try to set the 
    % fill value with NC_ATTPUT.
    if strcmp(attribute_name,'_FillValue')
        if library_version >= 4
            fmt = netcdf.inqFormat(ncid);
            switch(fmt)
                case {'FORMAT_CLASSIC','FORMAT_64BIT'}
                    % this is ok
                case {'FORMAT_NETCDF4','FORMAT_NETCDF4_CLASSIC'}
                    error('FillValues for netcdf-4 files should be set with NC_ADDVAR instead of NC_ATTPUT.');
            end

        end
    end

    % If netcdf-4 and the attribute is uint8, typecast it to int8.
    if isa(attval,'uint8');
        if library_version >= 4
            fmt = netcdf.inqFormat(ncid);
            switch(fmt)
                case {'FORMAT_CLASSIC','FORMAT_64BIT'}
                    % this is ok
                case {'FORMAT_NETCDF4','FORMAT_NETCDF4_CLASSIC'}
                    attval = typecast(attval,'int8');
            end

        end
    end

    if isnumeric(varname)
        varid = varname;
    else
        varid = netcdf.inqVarID(ncid, varname );
    end

    write_attribute(ncid,varid,attribute_name,attval);
    
    netcdf.endDef(ncid);

catch myException
    netcdf.close(ncid);
    rethrow(myException);
end

netcdf.close(ncid);

return;


%---------------------------------------------------------------------------
function write_attribute(ncid,varid,attribute_name,attval)

% If we are dealing with the fill value, then force the type to be
% correct.
if strcmp(attribute_name,'_FillValue')
    [name,xtype] = netcdf.inqVar(ncid,varid); %#ok<ASGLU>
    switch(xtype)
        case nc_double
            attval = double(attval);
        case nc_float
            attval = single(attval);
        case nc_int
            attval = int32(attval);
        case nc_short
            attval = int16(attval);
        case nc_byte
            attval = int8(attval);
        case nc_char
            attval = char(attval);
    end
    
    % If this is a netcdf-4 file, then don't treat _FillValue as an
    % attribute.  Let the library take care of it.
    v = netcdf.inqLibVers();
    if str2double(v(1)) > 3
        fmt = netcdf.inqFormat(ncid);
        if strcmp(fmt,'FORMAT_NETCDF4') || strcmp(fmt,'FORMAT_NETCDF4_CLASSIC')
            netcdf.defVarFill(ncid,varid,false,attval);
            return
        end
    end
end

try
    if iscellstr(attval) && (numel(attval) == 1)
        netcdf.putAtt(ncid,varid,attribute_name,attval{1});
    else
        netcdf.putAtt(ncid,varid,attribute_name,attval);
    end
catch me
    switch(me.identifier)
        case 'MATLAB:netcdf_common:emptySetArgument'
            % Bug #609383
            % Please consult the README.
            %
            % If char, change attval to ' '
            warning('SNCTOOLS:NCATTPUT:emptyAttributeBug', ...
                ['Changing attribute from empty to single space, ' ...
                'please consult the README regarding Bug #609383.']);
            netcdf.putAtt(ncid,varid,attribute_name,' ');
        otherwise
            rethrow(me);
    end
    
end
    
