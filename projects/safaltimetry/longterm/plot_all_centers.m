load('C:\Users\ashao\Data\allcenters.mat')

lon = 0:2:360;
window = 1;
ntrans = zeros(size(lon));
for i = length(lon)-1
    
    ntrans(i+1) = sum( centers(:,3)>=lon(i+1) & centers(:,3)<(lon(i+2)) );
    
end
    
    
    