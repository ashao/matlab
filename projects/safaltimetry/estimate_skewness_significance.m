inpath = '/ltraid4/aviso/alongtrack/sla/vxxc_matlab/';
outpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/bootstrap_skewness/';
files = dir([inpath 't*.mat'])
nfiles = length(files);
parpool(8)
%%

options = statset('UseParallel',true);
for tidx = 1:nfiles
    tic;
    fprintf('Track %d/%d...',tidx,nfiles)
    load([inpath files(tidx).name]);
    [ntime npts] = size(track.sla);
    track.skewnessci = zeros(npts,2);
    sla = track.sla;
    track.skewness = skewness(track.sla);
    calcpts = find(sum(~isnan(sla))>30);
    ncalc = length(calcpts);
    skewnessci = zeros(ncalc,2);
    fprintf('Number of Points: %d ...',ncalc)
    parfor ptidx = 1:ncalc
%         fprintf('Point %d/%d\n',ptidx,ncalc)
        data = sla(:,calcpts(ptidx));
        data = data(~isnan(data));
        skewnessci(ptidx,:) = bootci(100,{@skewness,data},'alpha',0.05);
        
    end
   
    
    track.skewnessci(calcpts,:) = skewnessci;
    track.sigidx = abs(sum(sign(track.skewnessci),2))==2;
    track = rmfield(track,{'sla','time'});
    save([outpath files(tidx).name],'track')
    fprintf('Completed in %e seconds\n',toc);
end