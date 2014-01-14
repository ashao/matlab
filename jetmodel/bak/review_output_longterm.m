function [  ] = review_output_longterm( inpath, outpath, fronts, start )

files=dir([inpath filesep '*.mat']);
nfiles=length(files);

for i=start:nfiles
    disp(['Track ' files(i).name])
    load([inpath filesep files(i).name],'array')
    lon_inter=[];
    lat_inter=[];
    notnan=find(~isnan(array.lon));
    for front=1:length(fronts)
        [templon templat]=intersections(fronts(front).lon,fronts(front).lat,array.lon(notnan),array.lat(notnan));
        lon_inter=[lon_inter ; templon];
        lat_inter=[lat_inter ; templat];
    end
    
    cont=true;
    while cont
        clf
        hold on
        centeropt=array.optimal(:,2);
        plot(array.lat,array.skewness)
        plot(centeropt,zeros(length(centeropt),1),'kx','LineWidth',2,'MarkerSize',10)
        plot(lat_inter,zeros(size(lat_inter)),'ko','Markersize',10,'LineWidth',2)
        xlim([-50 -20])
        [deletecenter null button]=ginput(1);
        if ~isempty(deletecenter) & button==1                        
            [null delidx]=min( (centeropt-deletecenter).^2);
            array.optimal(delidx,:)=[];
            save([outpath filesep files(i).name],'array')
            cont=true;
        else
            cont=false;
        end
        
    end    
end

