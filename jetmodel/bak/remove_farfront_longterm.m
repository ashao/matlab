function [  ] = remove_farfront_longterm( inpath, outpath, fronts, delta )


files=dir([inpath filesep '*.mat']);
nfronts=length(fronts);
npointsdel=0;


if ~exist(outpath,'dir')
    mkdir(outpath)
end

for i=1:length(files)

    load([inpath filesep files(i).name])
    [ncenters null]=size(array.optpar);
    delidx=[];
    notnan=find( ~isnan(array.lat));    
    dist=zeros(nfronts,1);
    
    for j=1:ncenters
        centerlat=array.optimal(j,2);
        centerlon=interp1(array.lat(notnan),array.lon(notnan),centerlat);
        for k=1:nfronts
            dist(k)=min( ((fronts(k).lat-centerlat)*111).^2 + ((fronts(k).lon-centerlon)*87).^2 );
        end
        dist=sqrt(min(dist));
        if dist>delta
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


end

