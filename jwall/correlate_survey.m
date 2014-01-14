importdata('jwall.survey.csv');
data = ans.data;
data(:,1)=[];
[nsamp nvars] = size(data);

[cvals pvals] = corrcoef(data);
% figure(1)
% pcolor(cvals)
% caxis([-1 1]);axis equal;xlim([1 nvars]);ylim([1 nvars]);colorbar;
% title('Correlation Coefficients');
% figure(2)
% pcolor(double(pvals<0.05))
% caxis([0 1]);axis equal;xlim([1 nvars]);ylim([1 nvars]);colorbar;
% title('Significant? 95% Level')

csvwrite('corrcoefficients.csv',cvals);
csvwrite('pvals.csv',pvals);
csvwrite('significance.csv',pvals<0.05)