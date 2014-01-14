function [ ] = auto_lag_endpoints( guesspath, inpath, outpath, start )

if nargin < 4
    start=1;
end
mkdir(outpath)
files=dir([guesspath filesep '*.mat']);
nfiles=length(files);

clf

latweight=0.8;
skewweight=-0.2;

for ii=start:nfiles
    
    disp(['Working on file ' inpath files(ii).name])
    load([guesspath filesep files(ii).name])
    if exist('track')
        array=track;
        clear track;
    end
    guessarray=array;
    
    load([inpath filesep files(ii).name])
    if exist('track')
        array=track;
        clear track;
    end            
    
    [ncenters null]=size(guessarray.optpar_pos);
    for jj=1:ncenters
                
        [array.optpar_pos(jj,1) array.optpar_pos(jj,2) array.optpar_pos(jj,3) ] = autoguess_center(array.lat,array.skewness_pos,guessarray.optpar_pos(jj,1),skewweight,latweight);
        [array.optpar_neg(jj,1) array.optpar_neg(jj,2) array.optpar_neg(jj,3) ] = autoguess_center(array.lat,array.skewness_neg,guessarray.optpar_neg(jj,1),skewweight,latweight);        
        
    end
   
    subplot(2,1,1); hold on;
    notnan=find(~isnan(array.skewness_pos));
    centerlon=interp1(array.lat(notnan),array.lon(notnan),array.optpar_pos(:,1));
    scatter(centerlon,array.optpar_pos(:,1))
    subplot(2,1,2); hold on;
    
    notnan=find(~isnan(array.skewness_pos));
    centerlon=interp1(array.lat(notnan),array.lon(notnan),array.optpar_neg(:,1));
    scatter(centerlon,array.optpar_neg(:,1))
    drawnow;
    save([outpath filesep files(ii).name],'array')
    
end

