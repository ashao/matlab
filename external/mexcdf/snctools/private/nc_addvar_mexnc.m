function nc_addvar_mexnc(ncfile,varstruct,preserve_fvd)


[ncid, status] = mexnc ( 'open', ncfile, nc_write_mode );
if ( status ~= 0 )
    ncerr = mexnc ( 'strerror', status );
    error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:OPEN', ...
        'OPEN failed on %s, ''%s''', ncfile, ncerr);
end

% determine the dimids of the named dimensions
num_dims = length(varstruct.Dimension);
dimids = zeros(num_dims,1);
for j = 1:num_dims
    [dimids(j), status] = mexnc ( 'dimid', ncid, varstruct.Dimension{j} );
    if ( status ~= 0 )
        mexnc ( 'close', ncid );
        ncerr = mexnc ( 'strerror', status );
        error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:DIMID', ncerr );
    end
end

% If preserving the fastest varying dimension in mexnc, we have to 
% reverse their order.
if preserve_fvd
    dimids = flipud(dimids);
end

status = mexnc ( 'redef', ncid );
if ( status ~= 0 )
    ncerr = mexnc ( 'strerror', status );
    mexnc ( 'close', ncid );
    error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:REDEF', ncerr );
end

% We prefer to use 'Datatype' instead of 'Nctype', but we'll try to be 
% backwards compatible.
if isfield(varstruct,'Datatype')
    [varid, status] = mexnc ( 'DEF_VAR', ncid, varstruct.Name, ...
        varstruct.Datatype, num_dims, dimids );
else
    [varid, status] = mexnc ( 'DEF_VAR', ncid, varstruct.Name, ...
        varstruct.Nctype, num_dims, dimids );
end
if ( status ~= 0 )
    ncerr = mexnc ( 'strerror', status );
    mexnc ( 'endef', ncid );
    mexnc ( 'close', ncid );
    error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:DEF_VAR', ncerr );
end


if ~isempty(varstruct.Chunking)

    if preserve_fvd
        chunking = fliplr(varstruct.Chunking);
    else
        chunking = varstruct.Chunking;
    end
    
    if ( numel(chunking) ~= num_dims) 
        mexnc ( 'endef', ncid );
        mexnc ( 'close', ncid );
        error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:defVarChunking', ...
           'Chunking size does not jive with number of dimensions.');
    end

    status = mexnc('DEF_VAR_CHUNKING',ncid,varid,'chunked',chunking);
    if ( status ~= 0 )
        ncerr = mexnc ( 'strerror', status );
        mexnc ( 'endef', ncid );
        mexnc ( 'close', ncid );
        error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:DEF_VAR_CHUNKING', ncerr );
    end
end

if (varstruct.Shuffle || varstruct.Deflate)

    status = mexnc('DEF_VAR_DEFLATE',ncid,varid, varstruct.Shuffle,varstruct.Deflate,varstruct.Deflate);
    if ( status ~= 0 )
        ncerr = mexnc ( 'strerror', status );
        mexnc ( 'endef', ncid );
        mexnc ( 'close', ncid );
        error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:DEF_VAR_DEFLATE', ncerr );
    end
end

status = mexnc ( 'enddef', ncid );
if ( status ~= 0 )
    ncerr = mexnc ( 'strerror', status );
    mexnc ( 'close', ncid );
    error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:ENDDEF', ncerr );
end

status = mexnc ( 'close', ncid );
if ( status ~= 0 )
    ncerr = mexnc ( 'strerror', status );
    error ( 'SNCTOOLS:NC_ADDVAR:MEXNC:CLOSE', ncerr );
end




% Now just use nc_attput to put in the attributes
for j = 1:length(varstruct.Attribute)
    attname = varstruct.Attribute(j).Name;
    attval = varstruct.Attribute(j).Value;
    nc_attput(ncfile,varstruct.Name,attname,attval);
end



%--------------------------------------------------------------------------
