function [retrieval_method,fmt] = snc_read_backend(ncfile)
%SNC_READ_BACKEND   determine which netCDF library to use
%
% which backend do we employ?  Many, many possibilities to consider here.
%
%   [retrieval_method,fmt] = snc_read_backend(ncfile)
%
% returns selection for specified file, http or url.
%
%   [retrieval_method,fmt] = snc_read_backend()
%
% returns all available retrieval_method and fmt options
%
%See also: snctools, snc_format

import ucar.nc2.dods.*    
import ucar.nc2.*
    
retrieval_methods.java     = 'java';
retrieval_methods.tmw_hdf4 = 'tmw_hdf4';
retrieval_methods.tmw_hdf4_2011a = 'tmw_hdf4_2011a';
retrieval_methods.mexnc    = 'mexnc';
retrieval_methods.tmw      = 'tmw';

fmts = snc_format();
fmt  = '';

if nargin==0
   retrieval_method = retrieval_methods;
else   
   retrieval_method = '';
end

% Check for this early.
if isa(ncfile,'ucar.nc2.NetcdfFile') 
    retrieval_method = retrieval_methods.java;
	fmt = fmts.netcdf_java;
	return
end

mv = version('-release');

fmt = snc_format(ncfile);

% These cases have no alternatives.
if strcmp(fmt,fmts.HDF4) 
	switch(mv)
		case {'14','2006a','2006b','2007a','2007b','2008a','2008b', ...
		      '2009a','2009b','2010a','2010b'} 
			  retrieval_method = retrieval_methods.tmw_hdf4;
		otherwise
			  retrieval_method = retrieval_methods.tmw_hdf4_2011a;
		end

    fmt = fmts.HDF4;
	return
elseif (strcmp(fmt,fmts.GRIB) || strcmp(fmt,fmts.GRIB2) || strcmp(fmt,fmts.URL))
    % Always use netcdf-java for grib files or URLs (when java is enabled).

    if ~exist('NetcdfFile','class')
        error('Netcdf-java must be available to read this.');
    end
    retrieval_method = retrieval_methods.java;
    fmt = fmts.netcdf_java;
	return
end

switch ( mv )
    case { '11', '12', '13' };
		error('Not supported on releases below R14.');

    case { '14', '2006a', '2006b', '2007a', '2007b', '2008a' }
		% No native matlab support here.  Use mexnc if available, 
		% otherwise try java.
        if strcmp(fmt,fmts.NetCDF)
			try
		    	v = mexnc('inq_libvers');
                retrieval_method = retrieval_methods.mexnc;
            catch
                if ~exist('NetcdfFile','class')
                    error('Either netcdf-java or the mexnc mex-file must be available in order to read this.');
                end
                retrieval_method = retrieval_methods.java;
            end
            fmt = fmts.NetCDF;
        elseif strcmp(fmt,fmts.NetCDF4)
            if ~exist('NetcdfFile','class')
                error('Netcdf-java must be available to read this.');
            end
            retrieval_method = retrieval_methods.java;
            % Last chance is if it is some format that netcdf-java can handle.
            % Hope for the best.
            fmt = fmts.NetCDF4;
        end
        
    case { '2008b', '2009a', '2009b', '2010a' }
        % 2008b introduced native netcdf-3 support.
        % netcdf-4 still requires either mexnc or java, and we will favor
        % java again.
        if strcmp(fmt,fmts.NetCDF)
            % Use TMW for all local netcdf-3 files.
            retrieval_method = retrieval_methods.tmw;
            fmt = fmts.NetCDF;
        elseif strcmp(fmt,fmts.NetCDF4)
            if ~exist('NetcdfFile','class') 
                error('Netcdf-java must be available to read this.');
            end
            retrieval_method = retrieval_methods.java;
            fmt = fmts.NetCDF4;
        else
            % not netcdf-3 or netcdf-4
            % Last chance is if it is some format that netcdf-java can handle.
            retrieval_method = retrieval_methods.java;
            fmt = fmts.netcdf_java;
        end

    otherwise
        % R2010b:  introduced netcdf-4 support.
        if strcmp(fmt,fmts.NetCDF) || strcmp(fmt,fmts.NetCDF4)
            retrieval_method = retrieval_methods.tmw;
            fmt = fmts.NetCDF;
        else
            % Last chance is if it is some format that netcdf-java can handle.
            if ~exist('NetcdfFile','class') 
                error('Netcdf-java must be available to read this.');
            end
            fmt = fmts.netcdf_java;
        end

end

if isempty(retrieval_method)
    error('SNCTOOLS:unknownBackendSituation', ...  
	      'Could not determine which backend to use with %s.  If the file format is not netCDF, the java backend must be enabled.', ...
       ncfile );
end
return
