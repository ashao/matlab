function [ ] = optcruncher_2011( inpath, outpath, tracknum )

files=dir([inpath filesep '*.mat']);
load([inpath filesep files(tracknum).name])
if exist('track','var')
    array=track;
    clear track
end

if ~exist(outpath,'dir')
    mkdir(outpath)
end

[ncenters null]=size(array.optpar);
outfile=[outpath filesep files(tracknum).name];    


if ~exist(outfile,'file')
	for ii=1:ncenters
	    center=array.optpar(ii,1);
	    minlat=array.optpar(ii,2);
	    maxlat=array.optpar(ii,3);
	    idx=find(array.lat>=minlat & array.lat<=maxlat);
	    simlat=array.lat(idx);
	    datskew=array.skewness(idx);     
    
	    if length(idx)<3 | center > -20 
	        array.optimal(ii,1:4)=NaN;
	    else
	        array.optimal(ii,:)=skewfit_ga(double(simlat),double(datskew),double(center));
	        width=sqrt(2*array.optimal(ii,1)+array.optimal(ii,3));
	        SNR=array.optimal(ii,4);
	        centeropt=array.optimal(ii,2);        
	        disp(sprintf('Track %d Center %d (Initial/Y0/Width/SNR): %f %f %f %f', ...
		tracknum,ii,center,centeropt,width,SNR))
	    end
	end
    
	save([outpath filesep files(tracknum).name],'array')
end
