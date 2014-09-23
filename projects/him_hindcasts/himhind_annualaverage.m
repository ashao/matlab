tempfile = '/home/ashao/uw-apl/data/offtrac/input/normalyear/temp.hind.nc';
saltfile = '/home/ashao/uw-apl/data/offtrac/input/normalyear/salt.hind.nc';
%%

salt = zeros(60,49,210,360);
temp = zeros(60,49,210,360);

sidx = 0;
for i = 1:60
    fprintf('Year %d/%d',i,60)
    start = [sidx 0 0 0];
    count = [12 inf inf inf];
    salt(i,:,:,:) = nanmean(nc_varget(saltfile,'salt',start,count));
%     temp(i,:,:,:) = nanmean(nc_varget(temp   file,'temp',start,count));
    sidx = sidx+12;
end