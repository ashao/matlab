inpath = '/ltraid4/weather_stations/zambia/';
files = dir([inpath '*.mat']);
nstations = length(files);
for i =1:nstations
   
    temp = load([inpath files(i).name]);
    fields = fieldnames(temp.station);
    
    for j = 1:length(fields)
       
        zambia(i).(fields{j}) = temp.station.(fields{j});
        
    end
    
end

temp = {zambia.prcp};

for i=1:nstations
   
    prcp(i) = temp{i};
    
end
%% Find overlapping time range

timerange(2) = inf;
timerange(1) = 0;
for i = 1:nstations
   
   timerange(1) = max([min(prcp(i).time) timerange(1)]);
   timerange(2) = min([max(prcp(i).time) timerange(2)]);
    
end

days = timerange(1):timerange(2);
precipitation.daily.data = nan([nstations,length(days)]);

for i = 1:nstations
    
    [null idxdays idxdata ] = intersect(days,prcp(i).time);
    precipitation.daily.data(i,idxdays) = prcp(i).data(idxdata);    
    
end
% 
% delidx = logical(sum(isnan(precipitation.daily.data)));
% days = days(~delidx);
precipitation.daily.time = days;
% precipitation.daily.data = precipitation.daily.data(:,~delidx);
%% Calculate monthly precipiation

[minyear null] = datevec(min(days));
[maxyear null] = datevec(max(days));
years = minyear:maxyear;
nyears = length(years);

precipitation.monthly.data = nan(nstations,nyears*12);
precipitation.monthly.time = nan(nstations,nyears*12);
[year month day] = datevec(precipitation.daily.time);
for i=1:nstations       
    idx = 1;
    for iyear = 1:nyears                        
        for mon = 1:12
            sumidx = year == years(iyear) & month == mon;
            precipitation.monthly.data(i,idx) = nanmean(precipitation.daily.data(i,sumidx));
            precipitation.monthly.time(i,idx) = datenum(years(iyear),mon,15);
            idx = idx + 1;
        end
    end
    
end
%% Calculate 3-month precipitation
precipitation.month3.data = nan(nstations,nyears*4);
precipitation.month3.time = nan(nstations,nyears*4);

for i = 1:nstations    
    sidx = 1;
    eidx = sidx+2;
    dataidx = 1;
    while eidx < 360
        precipitation.month3.data(i,dataidx) = nanmean(precipitation.monthly.data(i,sidx:eidx));
        precipitation.month3.time(i,dataidx) = nanmean(precipitation.monthly.time(i,sidx:eidx));
        dataidx = dataidx+1;
        sidx = eidx + 1;
        eidx = sidx+2;
    end    
end
% Calculate SPI for 3-month precip
precipitation.month3.spi = nan(nstations,nyears*4);
for i=1:nstations
   
    zeroidx = precipitation.month3.data == 0;
    notnan = ~isnan(precipitation.month3.data);
    data = precipitation.month3.data(~zeroidx & notnan);
    p = gamfit(data);
    incgam = gammainc(data./p(2),p(1));
    q = sum(zeroidx)/sum(notnan & ~zeroidx);
    h = q + (1-q)*incgam;
     
    precipitation.month3.spi(i,~zeroidx & notnan) = norminv(h,0,1);
end

%% Calculate 6-month precipitation
precipitation.month6.data = nan(nstations,nyears*2);
precipitation.month6.time = nan(nstations,nyears*2);

for i = 1:nstations    
    sidx = 1;
    eidx = sidx+5;
    dataidx = 1;
    while eidx <= 360
        precipitation.month6.data(i,dataidx) = nanmean(precipitation.monthly.data(i,sidx:eidx));
        precipitation.month6.time(i,dataidx) = nanmean(precipitation.monthly.time(i,sidx:eidx));
        dataidx = dataidx+1;
        sidx = eidx + 1;
        eidx = sidx+5;
    end    
end
% Calculate SPI for 6-month precip
precipitation.month6.spi = nan(nstations,nyears*2);
for i=1:nstations
   
    zeroidx = precipitation.month6.data(i,:) == 0;
    notnan = ~isnan(precipitation.month6.data(i,:));
    data = precipitation.month6.data(i,~zeroidx & notnan);
    p = gamfit(data);
    incgam = gammainc(data./p(2),p(1));
    q = sum(zeroidx)/sum(notnan & ~zeroidx);
    h = q + (1-q)*incgam;
     
    precipitation.month6.spi(i,~zeroidx & notnan) = norminv(h,0,1);
end
%% Calculate 12-month precipitation
precipitation.month6.data = nan(nstations,nyears);
precipitation.month6.time = nan(nstations,nyears);

for i = 1:nstations    
    sidx = 1;
    eidx = sidx+11;
    dataidx = 1;
    while eidx <= 360
        precipitation.month12.data(i,dataidx) = nanmean(precipitation.monthly.data(i,sidx:eidx));
        precipitation.month12.time(i,dataidx) = nanmean(precipitation.monthly.time(i,sidx:eidx));
        dataidx = dataidx+1;
        sidx = eidx + 1;
        eidx = sidx+11;
    end    
end
% Calculate SPI for 12-month precip
precipitation.month12.spi = nan(nstations,nyears);
for i=1:nstations
   
    zeroidx = precipitation.month12.data(i,:) == 0;
    notnan = ~isnan(precipitation.month12.data(i,:));
    data = precipitation.month12.data(i,~zeroidx & notnan);
    p = gamfit(data);    
    q = sum(zeroidx)/sum(notnan);
    h = q + (1-q)*gamcdf(precipitation.month12.data(i,notnan),p(1),p(2));
     
    precipitation.month12.spi(i,~zeroidx & notnan) = norminv(h,0,1);
end

corrmat = nan(nstations);
lagmat = nan(nstations);

for i=1:nstations
    for j=1:nstations
        [corrmat(i,j) lagmat(i,j)] = maxcorr( ...
            precipitation.month12.spi(i,:),precipitation.month12.spi(j,:),'coeff');
    end
end
corrmat(1:(nstations+1):end) = NaN;
lagmat(1:(nstations+1):end) = NaN;
