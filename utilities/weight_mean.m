function wmean = weight_mean(x,w);
% Calculates weighted mean of vector (x) with weights (w)
wmean=nansum(x.*w)./nansum(w);


end