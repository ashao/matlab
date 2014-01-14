function [out_pco2,out_year,hist_end] = f_co2_atmos(varargin)
% function [varargout]=f_co2_atmos(varargin)
%
% This function returns atmospheric values for CO2 for a certain year or a 
% range of years. The atmospheric conc. from the OCMIP-3/NOCES web site
% supplemented by recent data from the ESRL web-site are used to
% interpolate CO2 values to the years of the conc. Linear interpolation is used.
%
% Until 1990.0, the concentrations are based on spline fit to data from 
% ice cores (Siple) and Manua Loa (splmlo_co2_noces.dat). From 1990.5 to 2003.0 
% 12-month smoothed MLO (also Manua Loa?) are used (also splmlo_co2_noces.dat).
% After that concentrations from the annual mean Mauna Loa data reported at ESRL
% are supplemented (co2_annmean_mlo.txt).  (alternatively the latter dataset
% could started to be used in 1959 when the first annual mean is available). 
% 
% defaults used                         changeable by variable input argument
% -----------------------------------------------------------------------------
% range of years: 1765.5:0.5:end of sc. 'year' + <vector of years>
% years later than in tables set to NaN 'years_ahead' + <identifier> ('nan',
%                                       'constant','linear')
% second data set is used after 2003    'hist2_start' + <year>
%
% Output:
%
% out_pco2, out_year: column vectors of atmospheric pco2 and the corresponding
%                     year 
% hist_end: last year of historic data (=2010); from then on predictions 
%           depending on the scenario chosen may be used (predictions not
%           implemented yet)
%
% SM, 6/1/2011
%

% set defaults ================================================================
years_ahead = 'nan';
hist2_start_year = 2003.5;

% check variable input arguments and overwrite defaults ====================== 
i = 1;
while i <= length(varargin)

  if isstr(varargin{i})
    argin = upper(varargin{i});
    if strcmp(argin,'YEAR')
      out_year = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'YEARS_AHEAD')
      years_ahead = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'HIST2_START')
      hist2_start_year = varargin{i+1};
      i = i+2;
    else
      error_str=['f_co2_atmos: Unkown identifier ''' varargin{i} '''!'];
      error(error_str);
    end
  else
    error_str=['f_co2_atmos: Identifier must be a string!\n'];
    error(error_str);
  end              

end

% load the data with the atmospheric CO2 concentrations =======================
pco2_hist1 = load('/home/osprey/ashao/matlab/aco2/atmos_data/splmlo_co2_noces.dat');
yr1    = pco2_hist1(:,1);
pco2_1 = pco2_hist1(:,2);

pco2_hist2 = load('/home/osprey/ashao/matlab/aco2/atmos_data/co2_annmean_mlo_mod.txt'); % use modified file fitted for matlab ("%" at beginning of comment lines rather than "#")
yr2    = pco2_hist2(:,1)+0.5; % assume that annual mean represents mid-year
pco2_2 = pco2_hist2(:,2);
pco2_err_2 = pco2_hist2(:,3);

% merge dataset together ======================================================
i1 = find(yr1 < hist2_start_year - 1e-5);
i2 = find(yr2 >= hist2_start_year - 1e-5);
year = [yr1(i1); yr2(i2)];
pco2 = [pco2_1(i1); pco2_2(i2)];

% if no year column is given as input use years of data for output ============
if ~exist('out_year')
  out_year = year;
end

% find last year of historical data ===========================================
hist_end = max(year);

% take out duplicate entries for a year =======================================
% (this shouldn't actually happen the way the time histories are merged) 
%di = find(diff(year) == 0) + 1;
%year(di) = [];
%pco2(di) = [];

% interpolate atmospheric concentrations onto selected year ==================
% (data out of range gets set to NaN by the interpolation function) ==========

% initialize
out_pco2 = NaN*ones(size(out_year));

% interpolate data onto year vector  
gi = find(~isnan(out_year));
out_pco2(gi) = interp1(year,pco2,out_year(gi));

% set early atmoshperic concentrations before beginning of time history to =====
% const. value =================================================================
lyi = find(out_year < min(year));
out_pco2(lyi) = pco2(1) * ones(size(lyi));

% leave atmoshperic concentrations beyond years in table at NaN unless ========
% otherwise selected ==========================================================

switch upper(years_ahead)
case 'CONSTANT'
  [max_yr,maxi] = max(year);
  lyi = find(out_year > max_yr);
  out_pco2(lyi) = pco2(maxi) * ones(size(lyi));
case 'LINEAR'
  [max_yr,maxi] = max(year);
  lyi = find(out_year > max_yr);
  % not implemented yet
case 'NAN'
  % do not do anything because interpolation should have set values to NaN
end


