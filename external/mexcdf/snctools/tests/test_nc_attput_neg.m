function test_nc_attput_neg(ncfile,mode)

v = version('-release');
switch(v)
    case { '14','2006a','2006b','2007a'}
        fprintf('No negative tests run on %s...  ',v);
        return
    case {'2007b','2008a','2008b','2009a','2009b','2010a'}
		test_write_fill_value ( ncfile,mode );
        return
    otherwise
		test_write_fill_value ( ncfile,mode );
        return
end

return;







%--------------------------------------------------------------------------
function test_write_fill_value(ncfile,mode)
% Fill values for netcdf-4 files should only be set before the nc_enddef 
% has been called.  It's fine for netcdf-3 files, though.

ncfile = 'foo.nc';
nc_create_empty(ncfile,mode);
nc_adddim(ncfile,'x',5);
v.Name = 'y';
v.Dimension = {'x'};
nc_addvar(ncfile,v);

% If we try to write '_FillValue' as an attribute and if the file is netcdf-4,
% we should error out.
info = nc_info(ncfile);
if strcmp(info.Format,'NetCDF-4 Classic')

	% This should issue an error.
	try
		nc_attput(ncfile,'y','_FillValue',-99);
	catch
		return;
	end

	% If we get this far, then we failed to detect the condition.
   	error('failed');

else
	% OK it's not netcdf-4, so rewriting the fill value is fine.
	nc_attput(ncfile,'y','_FillValue',-99);
	fv = nc_attget(ncfile,'y','_FillValue');
	if fv ~= -99
    	error('failed');
	end

end
