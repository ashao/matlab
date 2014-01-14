function fix_creation_order_issue(ncfile)
% FIX_CREATION_ORDER_ISSUE
%   fix_creation_order_issue(nc4file) tries to fix certain netCDF-4 files 
%   by deleting the '_Netcdf4Dimid' attribute.  
%   


fid = H5F.open(ncfile,'H5F_ACC_RDWR','H5P_DEFAULT');

v = version('-release');
switch(v)
	case {'2011a', '2010b' };
		cleanup_group(fid,fid);
	case {'2010a','2009b'}
		cleanup_group_2010a(fid,'/');
	otherwise
	    cleanup_group(fid,fid);
end
H5F.close(fid);

try
	nc_dump(ncfile);
catch
	fprintf('No SNCTOOLS on your path, you need to manually verify that the issue is fixed.\n');
end


%--------------------------------------------------------------------------
function cleanup_group_2010a(parent_gid,child_group_name)

fprintf('Descending into group %s...\n', child_group_name);

child_gid = H5G.open(parent_gid,child_group_name);

idx = -1;
more_objects = true;
while more_objects
    idx = idx+1;
    objname = H5G.get_objname_by_idx(child_gid,idx);
    if isempty(objname)
        more_objects = false;
    else
        % get information about this group
        statbuf=H5G.get_objinfo(child_gid,objname,0);
        
        switch (statbuf.type)
            case H5ML.get_constant_value('H5G_GROUP')
                cleanup_group_2010a(child_gid,objname);
                
            case H5ML.get_constant_value('H5G_DATASET')
				obj_id = H5D.open(child_gid,objname);
                cleanup_object_2010a(obj_id);
				H5D.close(obj_id);
                
            case H5ML.get_constant_value('H5G_TYPE')
				obj_id = H5T.open(child_gid,objname);
                cleanup_object_2010a(obj_id);
				H5T.close(obj_id);
                
            otherwise
                fprintf ( '  Unhandled object type (does this ever happen? : %s\n', name);
        end
    end
end

H5G.close(child_gid);

fprintf('Ascending from group %s...\n', child_group_name);


%--------------------------------------------------------------------------
function cleanup_group(parent_gid,gid)

% netCDF uses indexing by creation order, so that's what we'll do too
idx_type = 'H5_INDEX_CRT_ORDER';
order = 'H5_ITER_INC';

child_group_name = H5I.get_name(gid);

fprintf('Descending into group %s...\n', child_group_name);

info = H5G.get_info(gid);

for j = 0:info.nlinks-1
	name = H5L.get_name_by_idx(parent_gid,child_group_name,idx_type,order,j,'H5P_DEFAULT');
    
	obj_id = H5O.open(gid,name,'H5P_DEFAULT');

    info = H5O.get_info(obj_id);
    switch(info.type)
    case H5ML.get_constant_value('H5O_TYPE_GROUP')
        cleanup_group(gid,obj_id);
    case H5ML.get_constant_value('H5O_TYPE_DATASET')
        cleanup_object(gid,name,obj_id);
    case H5ML.get_constant_value('H5O_TYPE_NAMED_DATATYPE')
        cleanup_object(gid,name,obj_id);
    end
    
    H5O.close(obj_id);

end

fprintf('Ascending from group %s...\n', H5I.get_name(gid));


%--------------------------------------------------------------------------
function cleanup_object(gid,objname,obj_id)

fprintf('Processing object %s... \n', objname);

% Does the attribute even exist?
try
    attrid = H5A.open_by_name(gid,objname,'_Netcdf4Dimid','H5P_DEFAULT','H5P_DEFAULT');
    % OK it's there.  Delete it.
    H5A.close(attrid);
    fprintf('\t''_Netcdf4Dimid'' removed from ''%s'' ...\n', objname);
    H5A.delete(obj_id,'_Netcdf4Dimid');
catch me
    % Most likely the attribute didn't exist.
    fprintf('\t''_Netcdf4Dimid'' does not exist, skipping ...\n');
end
    
return

%--------------------------------------------------------------------------
function cleanup_object_2010a(obj_id)
% 2010a lacks the "H5A.open_by_name" function

fprintf('Processing object %s... \n', H5I.get_name(obj_id));

% Does the attribute even exist?
try
    attrid = H5A.open_name(obj_id,'_Netcdf4Dimid');
    % OK it's there.  Delete it.
    H5A.close(attrid);
    fprintf('\t''_Netcdf4Dimid'' removed from ''%s'' ...\n', H5I.get_name(obj_id));
    H5A.delete(obj_id,'_Netcdf4Dimid');
catch me
    % Most likely the attribute didn't exist.
    fprintf('\t''_Netcdf4Dimid'' does not exist, skipping ...\n');
end
    
return
