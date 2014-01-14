%%
% WOA09 mapped fields for seasonally varying, nutrients, O2, etc are on 
% a vertical grid that that only goes down to 500m. So we use the annually
% averaged file which has 33 depth levels to fill in the rest.

indir = '/ltraid2/ashao/uw-apl/data/woa09/offtrac/';
files.annual = dir([indir filesep '*annual*.nc']);
files.monthly = dir([indir filesep '*monthly*.nc']);
field = {'o_an','n_an','p_an'};
nfields = length(field);
nmonths = 12; 

[nz nlat nlon] = size(nc_varget( ...
    [indir filesep files.annual(1).name],field{1}));
depth = nc_varget( ...
    [indir filesep files.annual(1).name],'depth');
time = nc_varget( ...
    [indir filesep files.annual(1).name],'time');

for iprop=1:nfields
    monarray = nc_varget( ...
        [indir filesep files.monthly(iprop).name],field{iprop});
    [ntime nz nlat nlon] = size(monarray);
    annarray(1,:,:,:) = nc_varget( ...
            [indir files.annual(iprop).name],field{iprop});
    for imon = 1:nmonths   
        fprintf('Property: %s Month %d\n',field{iprop},imon)        
        annarray(1,1:nz,:,:)=monarray(imon,:,:,:);        
        tempfile = [field{iprop} sprintf('.%02d.nc',imon)];
        copyfile(files.annual(iprop).name,tempfile);
        nc_varput(tempfile,field{iprop},annarray);
    end
end

%%
