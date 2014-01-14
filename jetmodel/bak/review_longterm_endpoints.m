function [  ] = review_longterm_endpoints( inpath, outpath, start )

files=dir([inpath filesep '*.mat']);
nfiles=length(files);

if ~exist(outpath,'dir')
    mkdir(outpath)
end

for i=start:nfiles
   
    load([inpath filesep files(i).name])
    if exist('track','var')
        array=track;
        clear track;
    end
    array
    [ncenters null]=size(array.optpar);
    
    for j=1:ncenters
        clf
        plot(array.lat,array.skewness)       
        hold on
        ylim([-1.5 1.5])
        xlim([-60 -40])
        center=array.optpar(j,1);
        minlat=array.optpar(j,2);
        maxlat=array.optpar(j,3);
        plot([minlat center maxlat],[1 0 -1])
        del=input('Delete Center?','s');
        if ~isempty(del)
            delidx=j;
        end
        
    end
    
    array.optpar(delidx,:)=[];
    save([outpath filesep files(i).name],'array');
    end
    
    
end


