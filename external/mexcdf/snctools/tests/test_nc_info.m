function test_nc_info(mode)

if nargin < 1
	mode = 'nc-3';
end

fprintf('\t\tTesting NC_INFO ...  ' );

run_negative_tests;

switch(mode)
	case 'nc-3'
		run_nc3_tests;
	case 'netcdf4-classic'
		run_nc4_tests;
	case 'netcdf4-enhanced'
		run_nc4_enhanced_tests;
	case 'http'
		run_http_tests;
	case 'opendap'
		run_opendap_tests;
end

fprintf('OK\n');

return


%--------------------------------------------------------------------------
function run_negative_tests()

test_no_inputs;
test_too_many_inputs;
test_file_not_netcdf;




%--------------------------------------------------------------------------
function run_nc4_enhanced_tests()

testroot = fileparts(mfilename('fullpath'));

ncfile = [testroot '/testdata/moons.nc'];
test_string_variable(ncfile);
test_global_string_attribute(ncfile);
test_empty_string_attribute(ncfile);


%--------------------------------------------------------------------------
function test_string_variable(ncfile)

info = nc_info(ncfile);
exp_data = struct('Name','ourano', ...
    'Nctype', 12, ...
    'Datatype', 'string', ...
    'Unlimited', 0, ...
    'Dimension', [], ...
    'Size', [2 3], ...
    'Attribute', [], ...
    'Chunking', [], ...
    'Shuffle', 0, ...
    'Deflate', 0 );
exp_data.Dimension = {'x' 'y'};
exp_data.Attribute = struct('Name','Bianca','Nctype',12,'Datatype','string', ...
    'Value',[]);
exp_data.Attribute.Value = {'Puck' 'Miranda'};

pfvd = getpref('SNCTOOLS','PRESERVE_FVD');
if pfvd
    exp_data.Dimension = {'y' 'x'};
    exp_data.Size = [3 2];
end
act_data = info.Dataset(1);
if ~isequal(act_data,exp_data)
    error('failed');
end

%--------------------------------------------------------------------------
function test_global_string_attribute(ncfile)

info = nc_info(ncfile);
act_data = info.Attribute;
exp_data = struct('Name','others', ...
    'Nctype', 12, ...
    'Datatype', 'string', ...
	'Value', []);
exp_data.Value = {'Francisco', 'Caliban', 'Stephano', 'Trinculo', ...
        'Sycorax', 'Margaret', 'Prospero', 'Setebos', 'Ferdinand'};
if ~isequal(act_data,exp_data)
    error('failed');
end

%--------------------------------------------------------------------------
function test_empty_string_attribute(ncfile)

info = nc_info(ncfile);
if ~isempty(info.Group.Dataset.Attribute(2).Value{1})
    error('failed');
end

%--------------------------------------------------------------------------
function run_nc4_tests()


testroot = fileparts(mfilename('fullpath'));

ncfile = [testroot '/testdata/empty-4.nc'];
test_emptyNetcdfFile(ncfile);

ncfile = [testroot '/testdata/just_one_dimension-4.nc'];
test_dimsButNoVars(ncfile);

ncfile = [testroot '/testdata/full-4.nc'];
test_smorgasborg(ncfile);

return




%--------------------------------------------------------------------------
function run_nc3_tests()


testroot = fileparts(mfilename('fullpath'));

ncfile = [testroot '/testdata/empty.nc'];
test_emptyNetcdfFile(ncfile);

ncfile = [testroot '/testdata/just_one_dimension.nc'];
test_dimsButNoVars(ncfile);

ncfile = [testroot '/testdata/full.nc'];
test_smorgasborg(ncfile);

return




%--------------------------------------------------------------------------
function test_no_inputs( )
try
	nc_info;
catch %#ok<CTCH>
    return
end
error ( 'succeeded when it should have failed.\n'  );





%--------------------------------------------------------------------------
function test_too_many_inputs()

testroot = fileparts(mfilename('fullpath'));
ncfile = fullfile(testroot, 'testdata/empty.nc');
try
	nc_info ( ncfile, 'blah' );
catch %#ok<CTCH>
    return
end
error('succeeded when it should have failed.');





%--------------------------------------------------------------------------
function test_file_not_netcdf()
ncfile = mfilename;
try
	nc_info ( ncfile );
catch %#ok<CTCH>
    return
end
error ( 'succeeded when it should have failed.' );







%--------------------------------------------------------------------------
function test_emptyNetcdfFile(ncfile)

nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	error( 'Filename was wrong.');
end
if ( ~isempty ( nc.Dimension ) )
	error( 'Dimension was wrong.');
end
if ( ~isempty ( nc.Dataset ) )
	error( 'Dataset was wrong.');
end
if ( ~isempty ( nc.Attribute ) )
	error('Attribute was wrong.');
end
return









%--------------------------------------------------------------------------
function test_dimsButNoVars(ncfile)

nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	error( 'Filename was wrong.');
end
if ( length ( nc.Dimension ) ~= 1 )
	error( 'Dimension was wrong.');
end
if ( ~isempty ( nc.Dataset ) )
	error( 'Dataset was wrong.');
end
if ( ~isempty ( nc.Attribute ) )
	error( 'Attribute was wrong.');
end
return










%--------------------------------------------------------------------------
function test_smorgasborg(ncfile)

nc = nc_info ( ncfile );
if ~strcmp ( nc.Filename, ncfile )
	error( 'Filename was wrong.');
end
if ( length ( nc.Dimension ) ~= 5 )
	error( 'Dimension was wrong.');
end
if ( length ( nc.Dataset ) ~= 6 )
	error( 'Dataset was wrong.');
end
if ( length ( nc.Attribute ) ~= 1 )
	error( 'Attribute was wrong.');
end
return






%--------------------------------------------------------------------------
function run_opendap_tests()

run_motherlode_test;
%--------------------------------------------------------------------------
function run_motherlode_test (  )

if getpref('SNCTOOLS','TEST_REMOTE',false) && ...
        getpref ( 'SNCTOOLS', 'TEST_OPENDAP', false ) 
    
    load('testdata/nc_info.mat');
    % use data of today as the server has a clean up policy
    today = datestr(floor(now),'yyyymmdd');
    url = ['http://motherlode.ucar.edu:8080/thredds/dodsC/satellite/CTP/SUPER-NATIONAL_1km/current/SUPER-NATIONAL_1km_CTP_',today,'_0000.gini'];
	fprintf('\t\tTesting remote DODS access %s...  ', url );
    
    info = nc_info(url);

    % Get rid of the time dependent attributes.
    info.Attribute(6:7) = [];
    
    if ~isequal(info,d.opendap.motherlode)
        error('failed');
    end
    fprintf('OK\n');
else
	fprintf('Not testing NC_DUMP on OPeNDAP URLs.  Read the README for details.\n');	
end
return





%--------------------------------------------------------------------------
function run_http_tests()

test_javaNcid;
return


%--------------------------------------------------------------------------
function test_javaNcid ()
import ucar.nc2.dods.*     
import ucar.nc2.*          

url = 'http://rocky.umeoce.maine.edu/GoMPOM/cdfs/gomoos.20070723.cdf';
jncid = NetcdfFile.open(url);
nc_info ( jncid );
close(jncid);
return



