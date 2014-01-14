function [ T,S,P, lat, lon ] = off_loadclim( datapath, month )
%OFF_LOADCLIM Load climatology fields for the given month
%   INPUT:
%       datapath (string) : path to offtrac input files
%       month: which month to extract will be mod(month,12), 1 corresponds
%           to January

if month>12
    monthdim=mod(month,12)-1;
    if monthdim==-1
        monthdim=11;
    end
end

tspath=[datapath filesep 'ts.nc'];
ppath=[datapath filesep 'gasx_ocmip2_himgrid.nc'];

tsdim=[1 4 210 360];
month_tab=[monthdim 0 0 0];
% monthdim
lon=nc_varget(tspath,'XH');
lat=nc_varget(tspath,'YH');
T=nc_varget(tspath,'TEMPCLIM',month_tab,tsdim);
S=nc_varget(tspath,'SALTCLIM',month_tab,tsdim);
P=nc_varget(ppath,'OCMIP_ATMP',[monthdim 0 0],[1 210 360]);

end

