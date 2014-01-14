function [ ] = fu_ryu_xsects(array,lon);

    lonidx=find(array.lonh==lon);
    
years=12*(5:5:60)-1;

for i=1:length(years)
    
    subplot(4,3,i)
    contourf(repmat(array.geolat(:,lonidx)',[30 1]), ...
        -squeeze(array.depth(years(i),:,:,lonidx)), ...
        double(squeeze(array.mn_cs137(years(i),:,:,lonidx))),20);
    xlim([5 40]);
    ylim([-600 0]);
    colorbar
    title(sprintf('137-Cs in Year %d along %3.1fW',2011+i*5,lon))
    caxis([0 5e-6])
    
end
    