function [centerlat minlat maxlat] = autoguess_center(lat, skew, guess, skewweight, latweight)

%% Truncate skewness record around guess and eliminate NaNs
idx=find(abs(lat-guess)<10);
lat=lat(idx);
skew=skew(idx);
arenan=find(isnan(skew));
lat(arenan)=[];
skew(arenan)=[];

%% Smooth the skewness curve
window=5;
filtskew=filtfilt(ones(window,1),1,skew)/window^2;

%% Loop until at least one candidate transition is found
delta=0.01;
ncans=0;
while ncans<1
    delta=delta-0.0005;
    
    %% Find peaks and valleys in the data
    [maxtab mintab]=peakdet(filtskew,delta,lat);    
        
    if ~isempty(maxtab)
        peaklats=maxtab(:,1);
        peakskew=maxtab(:,2);
        troughlats=mintab(:,1);
        troughskew=mintab(:,2);
        
        npeaks=length(maxtab);
        
        for ii=1:npeaks
            
            %% Find the trough associated with the skewness change
            troughidx=find( troughlats >  peaklats(ii) );
            
            if ~isempty(troughidx)
                ncans=ncans+1;
                [null, nexttroughidx]=min( troughlats(troughidx) );
                nexttroughidx=troughidx(nexttroughidx);
                
                latrange(ii,:)=[ peaklats(ii) troughlats(nexttroughidx) ];
                center(ii)=mean(latrange(ii,:));
                latdiff(ii)=abs(center(ii)-guess);
                skewdiff(ii)=abs(peakskew(ii)-troughskew(nexttroughidx));
            end
        end
    end
    
    
end
%% Scale latdiff and skewchange
latdiff=latdiff/max(latdiff);
skewdiff=skewdiff/max(skewdiff);

%% Evaluate penalty function for each potential transition
score=latweight*latdiff+skewweight*skewdiff;

[null, canidx]=min(score);

%% Output the transition candidate
centerlat=center(canidx);
minlat=latrange(canidx,1);
maxlat=latrange(canidx,2);

end
