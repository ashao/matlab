% function [data,prop,units] = f_add_alk0(data,prop,units,ocean,varargin)
%
% This function adds preformed alkalinity (ALK0) as a column following the
% empirical relationships given by Sabine et al. (2002) for the Pacific Ocean,
% Sabine et al. (1999) for the Indian Ocean, and Chung et al. (2003)  for 
% the Atlantic Ocean. If "global" is chosen as ocean, the only salinity-
% dependent relationship for Alk0 by Brewer et al. (1986) is used (as done
% for instance in Waugh et al. (2006) for estimating Cant). This can be changed 
% to use the formulas from the individual basins by passing "global_formula" + 
% "individual_ones" as variable input arguments (NOT IMPLEMENTED YET). In that
% case, the columns "   LON", "   LAT" must exist (changeable via "lon_col"/
% "lat_col" plus column name). 
%
% Also, the formula for the Atlantic Ocean may be varied to be the one by Lee et
% al. (2003) (very similar to Chung et al., 2003) by passing "atl_formula" + 
% "leeetal2003" as variable input arguments. "atl_formula","pac_formula", or 
% "ind_formula" + "breweretal1986" will change the formula for the individual
% basins to the Brewer et al. (1986) one (or 'breweretl1997' to Brewer et al.
% 1997)
% 
% The default names used for columns of salinity, oxygen, phosphate, nitrate, 
% and pot. temperature are 'SALNTY', 'OXYGEN', 'PHSPHT', 'NITRAT', and ' THETA'.
% The column added is named by default '  ALK0'. The dafault names can be 
% changed by passing vairbale input arguments. The identifiers 'alk0_col',
% 'sal_col','oxy_col','phspht_col', 'nitrat_col', 'theta_col' indicate that 
% that the next input argument is the name of the respective column that
% overwrites teh default. Depending on ocean only phosphate (Pacific, 
% Indian) or nitrate (Atlantic) are needed.  
%
% SM, 5/12/2011

function [data,prop,units] = f_add_alk0(data,prop,units,ocean,varargin)

% set default column descritpions =============================================

alk0_descr = '  ALK0';

sal_descr = 'SALNTY';
oxy_descr = 'OXYGEN';
phspht_descr = 'PHSPHT';
nitrat_descr = 'NITRAT';
theta_descr = ' THETA';

atl_formula = 'chungetal2003';
pac_formula = 'sabineetal2002';
ind_formula = 'sabineetal1999';
glob_formula = 'breweretal1986';

% check variable input arguments and replace defaults if argument passed ======
i = 1;
while i <= length(varargin)

  if isstr(varargin{i})
    argin = upper(varargin{i});
    if strcmp(argin,'ALK0_COL')
      alk0_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'SAL_COL')
      sal_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'OXY_COL')
      oxy_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'PHSPHT_COL')
      phspht_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'NITRAT_COL')
      nitrat_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'THETA_COL')
      theta_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'ATL_FORMULA') | strcmp(argin,'ATLANTIC_FORMULA')
      atl_formula = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'PAC_FORMULA') | strcmp(argin,'PACIFIC_FORMULA')
      pac_formula = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'IND_FORMULA') | strcmp(argin,'INDIAN_FORMULA')
      ind_formula = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'GLOB_FORMULA') | strcmp(argin,'GLOBAL_FORMULA') 
      glob_formula = varargin{i+1};
      i = i+2;
    else
      error_str=['f_add_alk0: Unknown identifier ''' varargin{i} '''!\n'];
      error(error_str);
    end
  else
    error_str=['f_add_alk0: Identifier must be a string!\n'];
    error(error_str);
  end

end

% find columns in data ========================================================

[m,n] = size(data);
csal = sel_col(prop,sal_descr);
coxy = sel_col(prop,oxy_descr);
cph  = sel_col(prop,phspht_descr);  % needed in Pac and Ind
cni  = sel_col(prop,nitrat_descr);  % needed in Atl 
cth  = sel_col(prop,theta_descr);

% leave this out since Brewer et al. (1986) formula only requires T & S
%if isempty(csal) |  isempty(coxy) | (isempty(cph) & isempty(cni)) | isempty(cth)
%  error_str = 'Not all necessary data columns available for calculating prefromed alkalintiy';
%  error(error_str);
%end

% calculate preformed alkalinity ==============================================
switch ocean
case {'atl','atlantic'}

  switch atl_formula
  case 'chungetal2003'
    % see Chung et al. (2003), GBC for formula used for Atlantic Ocean ========
    alk0 = 318.3 + 56.27*data(:,csal) + 0.09016*(data(:,coxy) + 10.625*data(:,cni));
  case 'leeetal2003'  % not surpringly very close to Chung et al. (2003)
    % see Lee et al. (2003), GBC for formula used for Atlantic Ocean ========
    alk0 = 335.7 + 55.80*data(:,csal) + 0.08924*(data(:,coxy) + 10.625*data(:,cni));
  case 'breweretal1986'  % greater than Chung (2003) for S>~33, less for S<33
    % see Brewer et al. (1986) as cited in Brewer et al. (1997) ==============
    alk0 = 50.56*data(:,csal) + 547.0;  % based on N. Atl. data
  case 'breweretal1997'  % lower than Chung (2003) *and* Brewer (1986) for 
                         % S>~33, greater for S<33 
    % see Brewer et al. (1997) ===============================================
    alk0 = 45.785*data(:,csal) + 703.7;  % based on N. Pac.  (WOCE P17) data
  end

case {'pac','pacific'}

  switch pac_formula
  case 'sabineetal2002'
    % see Sabine et al. (2002), GBC for formula used for Pacific Ocean ========
    alk0 = 148.7 + 61.36*data(:,csal) + 0.0941*(data(:,coxy) + 170*data(:,cph)) - 0.582*data(:,cth);
  case 'breweretal1986'
    % see Brewer et al. (1986) as cited in Brewer et al. (1997) ==============
    alk0 = 50.56*data(:,csal) + 547.0;  % based on N. Atl. data
  case 'breweretal1997'
    % see Brewer et al. (1997) ===============================================
    alk0 = 45.785*data(:,csal) + 703.7;  % based on N. Pac.  (WOCE P17) data
  end

case {'ind','indian'}

  switch ind_formula
  case 'sabineetal1999'
    % see Sabine et al. (1999), GBC for formula used for Indian Ocean ========
    alk0 = 378.1 + 55.22*data(:,csal) + 0.0716*(data(:,coxy) + 170*data(:,cph)) - 1.236*data(:,cth);
  case 'breweretal1986'
    % see Brewer et al. (1986) as cited in Brewer et al. (1997) ==============
    alk0 = 50.56*data(:,csal) + 547.0;  % based on N. Atl. data
  case 'breweretal1997'
    % see Brewer et al. (1997) ===============================================
    alk0 = 45.785*data(:,csal) + 703.7;  % based on N. Pac.  (WOCE P17) data
  end

case {'glob','global'}

  switch glob_formula
  case 'breweretal1986'
    % see Brewer et al. (1986) as cited in Brewer et al. (1997) ==============
    alk0 = 50.56*data(:,csal) + 547.0;  % based on N. Atl. data
  case 'breweretal1997'
    % see Brewer et al. (1997) ===============================================
    alk0 = 45.785*data(:,csal) + 703.7;  % based on N. Pac.  (WOCE P17) data
  case 'individual_ones'
    % not implemented yet
  end

end

data(:,n+1) = alk0;

% add property and units description of new column ============================

[mp,np] = size(prop);
pdescr = alk0_descr;
nn = length(pdescr);
if nn > np
  prop(mp+1,1:np) = pdescr(1:np);   % chop of rest of our description
else
  prop(mp+1,1:nn) = pdescr;
end

[mu,nu]=size(units);

udescr = 'umol/kg';

nn = length(udescr);

if nn > nu
  units(mu+1,1:nu) = udescr(1:nu);
else
  units(mu+1,1:nn) = udescr;
end

