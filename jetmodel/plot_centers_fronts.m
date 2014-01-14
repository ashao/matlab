[xypos xywidth]=postprocess_longterm(outpath,fronts([3 5]));
 xypos(:,3)=1:length(xypos);
 clf
 hold on
m_proj('Mercator','lon',[0 360],'lat',[-60 -30]);
m_plot(xypos(:,1),xypos(:,2),'x','MarkerSize',10,'LineWidth',2)
m_plot(fronts(3).lon,fronts(3).lat)
m_plot(fronts(5).lon,fronts(5).lat)
% m_text(xypos(:,1),xypos(:,2),num2str(xypos(:,3)))
m_coast('patch',[0 0 0]);
m_grid