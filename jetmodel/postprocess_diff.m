function [ posarray negarray diffarray ] = postprocess_diff( inpath )

files=dir([inpath filesep '*.mat']);
nl=length(files);
posarray=ones(5*nl,5)*NaN;
negarray=ones(5*nl,5)*NaN;

eidx=0;

for ii=1:nl
    
    load([inpath filesep files(ii).name])
    [ncenters(1) null]=size(array.optimal_pos);
    [ncenters(2) null]=size(array.optimal_neg);
    ncenters=min(ncenters);
    sidx=eidx+1;
    eidx=sidx+ncenters-1;
    notnan=find(~isnan(array.lat));
    posarray(sidx:eidx,1:4)=array.optimal_pos(1:ncenters,:);
    posarray(sidx:eidx,5)=interp1(array.lat(notnan),array.lon(notnan),array.optimal_pos(1:ncenters,2));
    negarray(sidx:eidx,1:4)=array.optimal_neg(1:ncenters,:);
    negarray(sidx:eidx,5)=interp1(array.lat(notnan),array.lon(notnan),array.optimal_neg(1:ncenters,2));
    
end
sidx=eidx+1;
posarray(sidx:end,:)=[];
negarray(sidx:end,:)=[];
posarray(:,6)=sqrt(2*posarray(:,1)+posarray(:,3).^2);
negarray(:,6)=sqrt(2*negarray(:,1)+negarray(:,3).^2);
% Rearrange so that lon/lat is first
posarray=posarray(:,[5 2 6 1 3 4]);
negarray=negarray(:,[5 2 6 1 3 4]);

diffarray=posarray-negarray;

bounds=0:60:360;
numshifts=find(abs(diffarray(:,2))>=0.2 & abs(diffarray(:,2))<=1);
disp(sprintf('Percentage of large shifts: %f',length(numshifts)/length(posarray)*100))
for ii=1:(length(bounds)-1);
    idx=find(posarray(:,1)>=bounds(ii) & posarray(:,1)<=bounds(ii+1));
    disp(sprintf('Longitude Range %d-%d mean shift: %f +/- %f',bounds(ii),bounds(ii+1),nanmean(diffarray(idx,2)),nanstd(diffarray(idx,2))/sqrt(length(idx))))
end



end
