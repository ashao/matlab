function station = parse_dly(filename)
%% station = parse_dly(filename)
% Parses the daily summaries of weather station data in the format used by
% GHCND and stores the variables, quality flags,
% Indices for the strings are based on
% ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt

column_ranges = {1:11,12:15,16:17,18:21};
for day = 1:31
    
    day_columns = cellfun(@(x) x+8*(day-1),{22:26,27,28,29},'UniformOutput',false);
    column_ranges = {column_ranges{:} day_columns{:}};
    
end

fid = fopen(filename);

%%

% dataline = fgetl(fid);
fid = fopen(filename);
datacounter = 1;
linenumber = 1;
dataline = fgetl(fid);
while ischar(dataline)
    
    station.id = dataline(column_ranges{1});
    year = str2double(dataline(column_ranges{2}));
    month = str2double(dataline(column_ranges{3}));
    property = lower(dataline(column_ranges{4}));
    
    colidx = 5;
    %     fprintf('Year %d Month %d Day %d\n',year,month,day);
    for day = 1:31
        
        
        
        data = str2double(dataline(column_ranges{colidx}));
        
        if data > -9999
            %         if isempty(mflag)
            station.(property).time(datacounter) = datenum(year,month,day);
            station.(property).data(datacounter) = data;
            colidx = colidx + 1;
            
            station.(property).mflag(datacounter) = dataline(column_ranges{colidx});
            colidx = colidx + 1;
            
            qflag = dataline(column_ranges{colidx});
            station.(property).qflag(datacounter) = qflag;
            colidx = colidx + 1;
            
            sflag = dataline(column_ranges{colidx});
            station.(property).sflag(datacounter) = sflag;
            colidx = colidx + 1;
            datacounter = datacounter + 1;
        end
        
        
    end
    
    dataline = fgetl(fid);
end

fclose(fid);
% end