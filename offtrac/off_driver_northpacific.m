% Driver program to extract and process North Pacific Offtrac output

%% Extract Offtrac Output
datapath='/scratch/data/offtrac/input/';
ncfile='/models/offtrac/gasex/clim_test.0863.nc';
months=1:852;
latrange=[25 45];
lonrange=[-200 -140];
layers=1:10;

array=off_procsec(datapath,ncfile,months,latrange,lonrange,layers);

%% Divide fields into yearly segments
[nmonths nlat nlon]=size(array.cfc11_relsat);
nyears=nmonths/12;
cfc11_sat_year=mat2cell(array.cfc11_relsat,repmat(12,[nyears 1]));
cfc12_sat_year=mat2cell(array.cfc12_relsat,repmat(12,[nyears 1]));
sf6_sat_year=mat2cell(array.sf6_relsat,repmat(12,[nyears 1]));

mldepth_year=mat2cell(array.mldepth,repmat(12,[nyears 1]));

%% Calculate minimum saturations and maximum depth
cfc11_minsat_year=cellfun(@min,cfc11_sat_year,'UniformOutput',false);
cfc12_minsat_year=cellfun(@min,cfc12_sat_year,'UniformOutput',false);
sf6_minsat_year=cellfun(@min,sf6_sat_year,'UniformOutput',false);
mldepth_max_year=cellfun(@max,mldepth_year,'UniformOutput',false);

cfc11_minsat_year=reshape(cell2mat(cfc11_minsat_year),[nyears nlat*nlon]);
cfc12_minsat_year=reshape(cell2mat(cfc12_minsat_year),[nyears nlat*nlon]);
sf6_minsat_year=reshape(cell2mat(sf6_minsat_year),[nyears nlat*nlon]);
mldepth_max_year=reshape(cell2mat(mldepth_max_year),[nyears nlat*nlon]);