function [ ] = plotfield( array, dt )

if nargin < 2
dt=[];
end

mask=load('../data/mask.mat');
mask=mask.mask;
[ntime numz numy numx]=size(array.F);

for j=1:10

for i=1:ntime
subplot(2,2,1)
contour( array.x, array.y, squeeze(array.salt(i,1,:,:)).*mask)
colorbar
title('Salinity')

subplot(2,2,2)
contour(array.x,array.y,squeeze(array.temp(i,1,:,:)).*mask)
colorbar
title('Temperature')

subplot(2,2,3)
contour(array.x,array.y,squeeze(array.F(i,1,:,:)).*mask)
colorbar
title('Solubility')

subplot(2,2,4)
contour(array.x,array.y,squeeze(array.flux(i,:,:)).*mask)
title('Gas Exchange Flux')
colorbar

pause(dt)

end

end
