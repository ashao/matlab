function [ lastidx ] = review_skewness_fit( inpath, outpath, start )

files=dir([inpath filesep '*.mat']);
npointsdel=0;


if ~exist(outpath,'dir')
    mkdir(outpath)
end

numtracks=1000;

for i=start:length(files)
    
    load([inpath filesep files(i).name])
    [ncenters null]=size(array.optimal);
    delidx=[];
    disp(sprintf('Track Number: %d',str2num(files(i).name(2:4))))    
    
    for j=1:ncenters
        
        minlat=array.optpar(j,2);
        maxlat=array.optpar(j,3);
        idx=find(array.lat>=minlat & array.lat<=maxlat);
        skewdat=array.skewness(idx);
        latdat=array.lat(idx);
        disp(sprintf('Original Center: %f',array.optpar(j,1)))
        disp(sprintf('Optimal Center: %f',array.optimal(j,2)))
        
        if ~isempty(latdat) & min(latdat)<-30
            randn_heights=make_zeroskewdist(numtracks,length(latdat));
            randn_x0=make_zeroskewdist(1,numtracks)';
            simskew=gensim(array.optimal(j,:),latdat,randn_heights,randn_x0);
            
            clf
            hold on
            plot(latdat,skewdat)
            plot(latdat,simskew,'LineWidth',2)
%             ylim([-1 1])
            grid on
            
            delchoice=input('Delete Center?','s');
            if ~isempty(delchoice)
                delidx=j;
            end
        else
            delidx=j;
        end
        
    end
    
    array.optimal(delidx,:)=[];
    array.optpar(delidx,:)=[];
    
    save([outpath filesep files(i).name]);
    
end

end

