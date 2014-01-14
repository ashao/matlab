function [ array ] = off_proc(datapath,ncfile,month,year,tracpath)
% OFF_PROC: Driver routine to extract SF6 and CFCs
% Calls off_extrday, off_conv2grav
% INPUT:
%   tspath (string): path to location of ts.nc
%   ncfile (string): path to Offtrac output
%   month: month to extract information
%   tracpath (string): path to tracer .mat files

% Load T,S,P
[T S P]=off_loadclim(datapath,month);

array.month=month;
array.lath=nc_varget([datapath 'metrics.nc'],'lath');
array.lonh=nc_varget([datapath 'metrics.nc'],'lonh');
array.geolat=nc_varget([datapath 'metrics.nc'],'geolat');
array.geolon=nc_varget([datapath 'metrics.nc'],'geolon');
array.wet=nc_varget([datapath 'metrics.nc'],'wet');
array.T=T;
array.S=S;
array.P=P;


% Extract concentration data from offtrac
array=off_extrday(ncfile,month,array);
array.sw_dens=sw_dens(S,T,array.depth);
array.sw_pden=sw_pden(S,T,array.depth,zeros(size(array.depth)));
thin_idx=find(array.mn_h<1e-3);
array.mn_cfc11(thin_idx)=NaN;
array.mn_cfc12(thin_idx)=NaN;
array.mn_sf6(thin_idx)=NaN;

% Convert to gravimetric units
% disp('Converting concentrations to gravimetric units');

array.cfc11_kg=real(off_conv2grav(array.mn_cfc11,T,S,array.depth));
array.cfc12_kg=real(off_conv2grav(array.mn_cfc12,T,S,array.depth));
array.sf6_kg=real(off_conv2grav(array.mn_sf6,T,S,array.depth));

% Calculate column inventory in gravimetric
% disp('Calculating column inventories')
array.cfc11_inv_kg=off_globalinv(array,array.cfc11_kg);
array.cfc12_inv_kg=off_globalinv(array,array.cfc12_kg);
array.sf6_inv_kg=off_globalinv(array,array.sf6_kg);

if nargin>3
    
    % Calculate relative saturation relative to surface saturation
%     disp('Calculating Relative Saturations')
    load([tracpath filesep 'cfc11.mat'])
    load([tracpath filesep 'cfc12.mat'])
    load([tracpath filesep  'sf6.mat'])
    
    month=mod(month,12);
    if month==0
        month=12;
    end
    
    ptT=sw_ptmp(S,T,array.depth,zeros(size(array.depth)));
    
    time=year+month/12;
    array.cfc11_sat=trac_calcsat( cfc11, time, ptT, S, P, array.lath, array.lonh );
    array.cfc12_sat=trac_calcsat( cfc12, time, ptT, S, P, array.lath, array.lonh );
    array.sf6_sat=trac_calcsat( sf6, time, ptT, S, P, array.lath, array.lonh );
    
    array.cfc11_relsat=array.mn_cfc11./array.cfc11_sat;
    array.cfc12_relsat=array.mn_cfc12./array.cfc12_sat;
    array.sf6_relsat=array.mn_sf6./array.sf6_sat;
    
    
end

end
