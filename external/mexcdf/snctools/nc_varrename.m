function nc_varrename ( ncfile, old_variable_name, new_variable_name )
%nc_varrename Rename netCDF variable.
%   nc_varrename(ncfile,oldvarname,newvarname) renames a netCDF variable
%   from OLDVARNAME to NEWVARNAME.
%
%   Example:
%      nc_create_empty('myfile.nc');
%      nc_adddim('myfile.nc','x',10);
%      v.Name = 'y';
%      v.Datatype = 'double';
%      v.Dimension = { 'x' };
%      nc_addvar('myfile.nc',v);
%      nc_dump('myfile.nc');
%      nc_varrename('myfile.nc','y','z');
%      nc_dump('myfile.nc');
%
%   See also nc_addvar.



backend = snc_write_backend(ncfile);
switch(backend)
	case 'tmw'
    	nc_varrename_tmw( ncfile, old_variable_name, new_variable_name )
	case 'mexnc'
    	nc_varrename_mexnc( ncfile, old_variable_name, new_variable_name )
    otherwise
    	error('NC_VARRENAME not supported with %s backend.', backend);
end



%--------------------------------------------------------------------------
function nc_varrename_hdf4(hfile,old_varname,new_varname )

fid = hdfh('open',hfile,'readwrite',0);
if fid < 0
    error('Could not open %s.', hfile);
end

status = hdfv('start',fid);
if status < 0
    hdfh('close',fid);
    error('Could not initialize Vgroup interface.');
end


% Look for the vgroup that represents the old variable.
vg_ref = hdfv('find',fid, old_varname);
if vg_ref < 0
    hdfv('end',fid);
    hdfh('close',fid);
    error('Could not find SDS %s.', old_varname);
end

% Get access to the vgroup 
vg_id = hdfv('attach',fid, vg_ref, 'w');
if vg_id < 0
    hdfv('end',fid);
    hdfh('close',fid);
    error('Could not attach to Vgroup for %s.', old_varname);
end

% Change from "SDS A" to "SDS B" 
status = hdfv('setname', vg_id, new_varname);
if status < 0
    hdfv('detach',vg_id);
    hdfv('end',fid);
    hdfh('close',fid);
    error('Could not rename variable from %s to %s.', old_varname, new_varname);
end

% Terminate access to the vgroup. 
status = hdfv('detach',vg_id);
if status < 0
    hdfv('end',fid);
    hdfh('close',fid);
    error('Could not detach from vgroup.');
end


status = hdfv('end',fid);
if status < 0
    hdfh('close',fid);
    error('Could not close down Vgroup interface.');
end

hdfh('close',fid);

%--------------------------------------------------------------------------
function nc_varrename_mexnc ( ncfile, old_variable_name, new_variable_name )
[ncid,status ]=mexnc('OPEN',ncfile,nc_write_mode);
if status ~= 0
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARGET:MEXNC:OPEN', ncerr );
end


status = mexnc('REDEF', ncid);
if status ~= 0
    mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARGET:MEXNC:REDEF', ncerr );
end


[varid, status] = mexnc('INQ_VARID', ncid, old_variable_name);
if status ~= 0
    mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARGET:MEXNC:INQ_VARID', ncerr );
end


status = mexnc('RENAME_VAR', ncid, varid, new_variable_name);
if status ~= 0
    mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARGET:MEXNC:RENAME_VAR', ncerr );
end


status = mexnc('ENDDEF', ncid);
if status ~= 0
    mexnc('close',ncid);
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARGET:MEXNC:ENDDEF', ncerr );
end


status = mexnc('close',ncid);
if status ~= 0
    ncerr = mexnc('strerror', status);
    error ( 'SNCTOOLS:NC_VARGET:MEXNC:CLOSE', ncerr );
end


