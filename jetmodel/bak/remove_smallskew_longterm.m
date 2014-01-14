function [  ] = remove_smallskew_longterm( inpath, outpath, delta )

files=dir([inpath filesep '*.mat']);
npointsdel=0;

if ~exist(outpath,'dir')
    mkdir(outpath)
end

for i=1:length(files)

    load([inpath filesep files(i).name])
    [ncenters null]=size(array.optpar);
    delidx=[];
    for j=1:ncenters
        minlat=array.optpar(j,2);
        maxlat=array.optpar(j,3);
        idx=find(array.lat >= minlat & array.lat <= maxlat);
        diffskew=max(array.skewness(idx))-min(array.skewness(idx));
        if diffskew < delta
           delidx=[delidx j];           
        end
    end
    
    array.optpar(delidx,:)=[];
    array.optimal(delidx,:)=[];
    npointsdel=npointsdel+length(delidx);
%     disp(sprintf('Track %s deleted %d points',files(i).name,length(delidx)))
    
    save([outpath filesep files(i).name],'array')
end

disp(sprintf('Deleted %d points',npointsdel))