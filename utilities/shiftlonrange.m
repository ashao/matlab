function [ shiftlon shiftdata ] = shiftlonrange(lon,data,maxlon)
[m n ] = size(lon);
if m==1
   lon = lon'; 
end

shiftcol = min(find(lon>maxlon))-1;
% size(lon)
shiftlon = circshift(lon,-shiftcol);
shiftlon(shiftlon>maxlon)=shiftlon(shiftlon>maxlon)-360;
shiftdata= circshift(data,[0 -shiftcol]);

end