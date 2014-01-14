function [minlat maxlat]=set_endpoints(lat,skew,centerlat)


    arenan=find(isnan(skew));
    skew(arenan)=[];
    lat(arenan)=[];

    filtskew=filtfilt(ones(3,1),1,skew)/3^2;
%    [m n]=size(track.optpar)r
%    for i=1:m
%	centerlat=track.optpar(i,2);
	            [maxtab mintab]=peakdet(filtskew,0.001,lat);

            peak_lat=maxtab(:,1);
            valley_lat=mintab(:,1);

            peak_candidate=peak_lat(find(peak_lat<=centerlat));
            valley_candidate=valley_lat(find(valley_lat>=centerlat));

            minlat=max(peak_candidate);
            maxlat=min(valley_candidate);
%	end

end

