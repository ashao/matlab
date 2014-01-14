function fileinfo = nc_info_tmw ( ncfile )
% NC_INFO Backend for Mathworks package.

ncid=netcdf.open(ncfile, nc_nowrite_mode );

fileinfo = nc_group_info_tmw(ncfile,ncid);

v = netcdf.inqLibVers;
if v(1) == '4'
    switch(netcdf.inqFormat(ncid))
        case 'FORMAT_CLASSIC'
            fileinfo.Format = 'NetCDF-3 Classic';
        case 'FORMAT_64BIT'
            fileinfo.Format = 'NetCDF-3 64bit';
        case 'FORMAT_NETCDF4_CLASSIC'
            fileinfo.Format = 'NetCDF-4 Classic';
        case 'FORMAT_NETCDF4'
            fileinfo.Format = 'NetCDF-4';
    end
end

netcdf.close(ncid);


fileinfo.Filename = ncfile;

%--------------------------------------------------------------------------
function info = nc_group_info_tmw(ncfile,ncid)

[ndims, nvars, ngatts] = netcdf.inq(ncid);

v = netcdf.inqLibVers;
if v(1) == '3'
    dimids = 0:ndims-1;
    info.Name = '/';
else
    dimids = netcdf.inqDimIDs(ncid);
    info.Name = netcdf.inqGrpNameFull(ncid);
end

% Get the dimensions
if ndims == 0
	Dimension = struct ( [] );
else
    Dimension = struct('Name','','Length',[],'Unlimited',false);
    Dimension = repmat(Dimension,ndims,1);
	for j = 1:ndims
		Dimension(j)=nc_getdiminfo_tmw(ncid,dimids(j));
	end
end



% Get the global attributes.
if ngatts == 0
	info.Attribute = struct([]);
else
    Attribute = struct('Name','','Nctype','','Datatype','','Value',NaN);
	for attnum = 0:ngatts-1
		Attribute(attnum+1) = nc_getattsinfo_tmw(ncfile,ncid,nc_global,attnum);
	end
	info.Attribute = Attribute;
end


% Get the variable information.
if nvars == 0
	Dataset = struct([]);
else
    Attribute = struct('Name','','Nctype','','Datatype','','Value',NaN);
    Dataset = struct('Name','','Nctype','','Datatype','','Unlimited',false,'Dimension',{''},'Size',[],'Attribute',Attribute,'Chunking',[],'Shuffle',[],'Deflate',[]);
	Dataset = repmat ( Dataset, nvars, 1 );
	for varid=0:nvars-1
		Dataset(varid+1) = nc_getvarinfo_tmw(ncfile,ncid,varid);
	end
end

info.Dimension = Dimension;
info.Dataset = Dataset;

Group = [];
if v(1) == '4'
    % Any groups?
    fmt = netcdf.inqFormat(ncid);
    if strcmp(fmt,'FORMAT_NETCDF4')
        childGroups = netcdf.inqGrps(ncid);
        if numel(childGroups) > 0
            Group = nc_group_info_tmw(ncfile,childGroups(1));
            Group = repmat(Group, numel(childGroups),1);
            for j = 2:numel(childGroups)
                Group(j) = nc_group_info_tmw(ncfile,childGroups(j));
            end
        end
    end
end
info.Group = Group;


return
