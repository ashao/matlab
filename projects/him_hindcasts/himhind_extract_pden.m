inpath = '/ltraid3/ashao/uw-apl/models/HIM/hindcast/';
files = dir([inpath 'sigmatheta*.mat']);
nyears = length(files);
him.pden = zeros(nyears*12,210,360),'single';
sidx = 1;
for i = 1:nyears
    
    fprintf('Year %d/%d\n',i,nyears);
    load([inpath files(i).name]);
    eidx = sidx+11;
    him.pden(sidx:eidx,:,:) = single(pden(:,1,:,:));
    sidx = eidx+1;
    
end