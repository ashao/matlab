paths.core = '/ltraid4/ashao/COREv2/data_IAF/CORRECTED/combined_years/';
paths.hind = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/';
outpath = '/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/';
infiles.u10 = [paths.core 'u_10.1948-2007.06JUN2011.nc'];
infiles.v10 = [paths.core 'v_10.1948-2007.06JUN2011.nc'];
infiles.slp = [paths.core 'slp.1948-2007.06JUN2011.nc'];
syear = 1948;
eyear = 2007;
years = syear:eyear;
nyears = length(years);

nmonths = nyears*12;
fillmat3d = zeros(nmonths,210,360);
fillmat4d = zeros(nmonths,49,210,360);

variables4d = {'temp','salt'};
variables3d = {'p_surf','CN','wind'};

eomdays = repmat(eomday(2002,1:12),[1 nyears]);
time = eomdays./2 + cumsum(eomdays)-31;

%% Create the output file 
% delete(outfile)
load metrics

dim3d = {'xh',360,'yh',210,'Time',nmonths};
dim4d = {'xh',360,'yh',210,'zl',49,'Time',nmonths};

% Do the 4D variables first
for i=1:length(variables4d)
    
    name = variables4d{i};
    outfile = [outpath name '.hind.nc'];
    if exist(outfile,'file')
        delete(outfile)
    end
    nccreate(outfile,'Time','Dimensions',{'Time',nmonths},'Format','64bit');
    nccreate(outfile,'zl','Dimensions',{'zl',49});
    nccreate(outfile,'xh','Dimensions',{'xh',360});
    nccreate(outfile,'yh','Dimensions',{'yh',210});
    ncwrite(outfile,'Time',time);
    ncwrite(outfile,'zl',1:49);
    ncwrite(outfile,'xh',metrics.lonh.data);
    ncwrite(outfile,'yh',metrics.lath.data);
    ncwriteatt(outfile,'/','creation_date',datestr(now));    
    fprintf('Adding %s...\n',name);
    nccreate(outfile,name,'Dimensions',dim4d);
end

% Now for the 3D variables
for i=1:length(variables3d)
    name = variables3d{i};
    outfile = [outpath name '.hind.nc'];
    if exist(outfile,'file')
        delete(outfile)
    end
    nccreate(outfile,'Time','Dimensions',{'Time',nmonths},'Format','64bit');    
    nccreate(outfile,'xh','Dimensions',{'xh',360});
    nccreate(outfile,'yh','Dimensions',{'yh',210});
    ncwrite(outfile,'Time',time);    
    ncwrite(outfile,'xh',metrics.lonh.data);
    ncwrite(outfile,'yh',metrics.lath.data);
    ncwriteatt(outfile,'/','creation_date',datestr(now));    
    fprintf('Adding %s...\n',name);
    nccreate(outfile,name,'Dimensions',dim3d);
end
%% Start adding in the data

% Do the easy stuff first 

outfiles.temp = [outpath 'temp.hind.nc'];
outfiles.salt = [outpath 'salt.hind.nc'];
outfiles.ice = [outpath 'CN.hind.nc'];

for i=1:nyears
    
    year = years(i);
    fprintf('Writing year %d...',year);
    start4d=[(i-1)*12  0 0 0 ];    
    start3d=[ (i-1)*12 0 0  ];    
    count3d = [12 210 360];
    count4d = [12 49 210 360];
    
    infiles.ocean = [paths.hind filesep sprintf('ocean_month.%d.nc',year)];
    infiles.ice = [paths.hind filesep sprintf('ice_month.%d.nc',year)];
    mondata.temp = nc_varget(infiles.ocean,'temp');
    mondata.salt = nc_varget(infiles.ocean,'salt');
    mondata.ice = squeeze(sum(nc_varget(infiles.ice,'CN'),2));       
%     
%     mondata.temp = permute(mondata.temp,[ 4 3 2 1 ]);
%     mondata.salt = permute(mondata.salt,[ 4 3 2 1 ]);
%     mondata.ice = permute(mondata.ice,[ 3 2 1 ]);
    fprintf('temp...');
    nc_varput(outfiles.temp,'temp',mondata.temp,start4d,count4d);
    fprintf('salt...');
    nc_varput(outfiles.salt,'salt',mondata.salt,start4d,count4d);
    fprintf('ice...');
    nc_varput(outfiles.ice,'CN',mondata.ice,start3d,count3d);   
    fprintf('Done!\n')
    
end

%% Now do the stuff that needs averaging and regridding...boo :(
load metrics.mat
clear fillmat4d mondata
eomdays = cumsum(repmat(eomday(2002,1:12),[1 nyears]));
somdays = cumsum(repmat(eomday(2002,1:12),[1 nyears]))-31;
time = ncread(infiles.slp,'TIME');
templon = nc_varget(infiles.slp,'LON');
templat = nc_varget(infiles.slp,'LAT');
outfiles.wind = [outpath 'wind.hind.nc'];
outfiles.p_surf = [outpath 'p_surf.hind.nc'];

fillmat3d = zeros(12,210,360);

mondata.wind = fillmat3d;
mondata.slp = fillmat3d;


templat(templat==max(templat(:)))=max(metrics.geolat.data(:));
[templon templat] = meshgrid(templon,templat);

templon = [fliplr(templon-360) templon];
templat = [fliplr(templat) templat];

yeartime = 0:0.25:365;
somdays = [0 cumsum(eomday(2001,1:11))];
eomdays = cumsum(eomday(2001,1:12));

for i=1:nyears
   fprintf('Writing year %d Month: ',years(i));
   % Grab the entire year's data
   shifttime = (i-1)*365;
   stime = 0 + shifttime;
   etime = stime + 366;   
   
   yearidx = find(time>stime & time<etime);
   start3d = [min(yearidx)-1 0 0];
   count3d = [length(yearidx) inf inf ];
   yeardata.u10 = nc_varget(infiles.u10,'U_10_MOD',start3d,count3d);
   yeardata.v10 = nc_varget(infiles.v10,'V_10_MOD',start3d,count3d);
   yeardata.slp = nc_varget(infiles.slp,'SLP',start3d,count3d);
   
   for mon = 1:12
       
       fprintf('%d...',mon);       
       avgidx = find(yeartime>=somdays(mon) & yeartime<eomdays(mon));
       wind_temp = sqrt(yeardata.u10.^2+yeardata.v10.^2);       
       wind_temp = squeeze(mean(wind_temp(avgidx,:,:)));
       wind_temp = [fliplr(wind_temp) wind_temp];
       mondata.wind(mon,:,:) = griddata(templon,templat,wind_temp, ...
           metrics.geolon.data, metrics.geolat.data);              
       
       p_temp = squeeze(mean(yeardata.slp(avgidx,:,:)));
       p_temp = [fliplr(p_temp) p_temp];
       mondata.slp(mon,:,:) = griddata(templon,templat,p_temp, ...
           metrics.geolon.data, metrics.geolat.data);
   end
   
   fprintf('DONE!\n')
   start3d = [ (i-1)*12 0 0];
   count3d = [12 210 360];
   nc_varput(outfiles.wind,'wind',mondata.wind,start3d,count3d);
   nc_varput(outfiles.p_surf,'p_surf',mondata.slp,start3d,count3d);
   
end