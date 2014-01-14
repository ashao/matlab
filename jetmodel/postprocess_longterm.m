function [ xypos xywidth ] = postprocess_longterm( inpath, fronts )

    files=dir([inpath filesep '*.mat']);
    nfiles=length(files);
    
    outarray=zeros(10*nfiles,6)*NaN;
    
    sidx=1;
    for ii=1:nfiles
        
        load([inpath filesep files(ii).name])
        [ncenters null]=size(array.optimal);
        eidx=sidx+ncenters-1;
        outarray(sidx:eidx,1:4)=array.optimal;
        notnan=find(~isnan(array.lat) & ~isnan(array.lon));
        outarray(sidx:eidx,5)=interp1(array.lat(notnan),array.lon(notnan),array.optimal(:,2));
        outarray(sidx:eidx,6)=str2num(files(ii).name(2:4));
        sidx=eidx+1;

    end
    
    arenan=find( isnan(outarray(:,1)) | isnan(outarray(:,5)) );
    outarray(arenan,:)=[];
    xypos=[outarray(:,5) outarray(:,2)];
    xywidth(:,1:2)=xypos;
    xywidth(:,3)=sqrt(outarray(:,1)*2+outarray(:,3).^2);
    xywidth(:,4)=outarray(:,6);

end

