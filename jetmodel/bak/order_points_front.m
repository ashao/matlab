function [ ordertab orderidx ] = orderpoints( xytab,fronts,startidx )

lon=0;
lat=0;
button=0;
ordertab=[];
count=startidx;
prevlon=0;

m_proj('Mercator','lon',[0 360],'lat',[-60 -30]);
subplot(2,1,2)
hold on
m_coast('patch',[0 0 0]);
m_grid
m_plot(xytab(:,1),xytab(:,2),'o')

[ncenters nparams]=size(xytab)

while ~isempty(lon+lat) & button ~=3    
    
    
    subplot(2,1,1,'replace')
    hold on
    for frontnum=1:length(fronts)
        plot(fronts(frontnum).lon,fronts(frontnum).lat)
    end        
    
    if nparams>2
        text(xytab(:,1),xytab(:,2),num2str(xytab(:,3),2));
    end
    
    plot(xytab(:,1),xytab(:,2),'x','LineWidth',2,'MarkerSize',10)
    if ~isempty(ordertab)
        plot(ordertab(:,1),ordertab(:,2))
    end
    xlim([prevlon-20 prevlon+20])
    [lon lat button]=ginput(1);
%     [lon lat]=m_xy2ll(xproj,yproj);    
    dist=zeros(length(xytab),1)*Inf;
    
    
    for j=1:length(xytab)
        dist(j)=m_lldist([lon xytab(j,1)],[lat xytab(j,2)]);
    end
    
    [null minidx]=min(dist);
    orderidx(count)=minidx;
    ordertab(count,:)=xytab(minidx,:);
    prevlon=ordertab(count,1);
    subplot(2,1,2)
    m_plot(ordertab(count,1),ordertab(count,2),'x','LineWidth',2,'MarkerSize',10);
    
    count=count+1;
    
    
    
end

end