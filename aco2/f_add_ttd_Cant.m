% function [data,prop,units] = f_add_ttd_Cant(data,prop,units,varargin)
%
% This function adds anthropogenic CO2 concentrations as a new column using 
% an inverse Gaussian TTD following Waugh et at. (2006). Columns for mean age
% (Gamma) and width (Delta) of the TTD must exist in the data matrix.
%
% input:   data   = data matrix (mxn) 
%          prop   = property martix (nx6) describing the data columns  
%          units  = units martix (nx7) describing the units of the data columns
%
% defaults used                       changeable by variable input arguments
% -----------------------------------------------------------------------------
% 'T12AGE' for mean TTD age column       'ttd_mean_age_col', <new name>
% 'T12WID' for width of TTD column       'ttd_width_col',,new name>      
% 'SALNTY' for salinity column           'sal_col', <new name>
% ' THETA' for temperature column        'temp_col', <new name>
% '  ALK0' for preformed alkalinity col. 'alk0_col', <new name>
% '  DATE' for date col.                 'date_col', <new name>
%
% '  CANT' for new Cant column           'cant_col', <column name>
%
% Error columns to propagate uncertainties  (NOT IMPLEMENTED YET):
% 'AGEERR' for TTD age error column        'age_err_col', <column name>
% 'AGEERL' for TTD age low error column    'age_erl_col', <column name>
% 'AGEERH' for TTD age high error column   'age_erh_col', <column name>
% 'WIDERR' for TTD width error column      'width_err_col', <column name>
% 'WIDERL' for TTD width low error column  'width_erl_col', <column name>
% 'WIDERH' for TTD width high error column 'width_erh_col', <column name>
% 'CANTER' for new Cant error column       'cant_err_col', <column name>
%
% atmos. CO2 concentrations = splmlo_co2_noces.dat + mlo annual mean from ESRL
%                                      'atmos_co2',<name of CO2 history> 
%                                       NOT IMPLEMENTED YET
%  surface satuturation = 100%         'surface_sat',<Co2 surface saturation> 
%                                       (= scale fatcor for atmospheric CO2)
%                                        NOT IMPLEMENTED YET
%  include_error = 0                   'include_error',<flag> NOT IMPLEMENTED YET
%  meas_error_on = 1                   'meas_error_on',<flag> NOT IMPLEMENTED YET
%  sol_error_on = 1                    'sol_error_on',<flag> NOT IMPLEMENTED YET
%  sat_error_on = 0                    'sat_error_on',<flag> NOT IMPLEMENTED YET
%  atmos_error_on = 1                  'atmos_error_on',<flag>; this assumes
%                                      that files exist with _errl and _errh
%                                      extensions run with high and low
%                                      atmos values  NOT IMPELMENTED YET
%
% '  DATE' for column with time of measurement
%  (in year-fraction format)           'date_col', <column name>
%
%  ----                                'meas_year', time of measurement in
%                                      year-fraction format: only one value for
%                                     whole cruise which will be used instead
%                                     of looking for a date column  
%                                     IMPLEMENTED YET
% hemisphere = 'N'                     'hemisphere',<indicator which hemisphere>
%                                      ('N', 'S', 'NS') NOT IMPLEMENTED YET  
%  '   LAT' for latitude column       'lat_col',<new name>; only used to identify
%                                     hemisphere if hemisphere = 'NS' chosen; NOT
%                                     IMPLEMENTED YET
% interp_lat_range = []            'interp_lat_range',2x1 vector of latitudes
%                                     for which interpolation between SH&NH
%                                     should happen; only used to identify
%                                     hemisphere if hemisphere = 'NS' chosen
%                                    NOT IMPELMENTED YET
%
% output:  data   = data matrix (mxn+2)       |
%          prop   = property martix (n+2x6)   | - 2 new columns added
%          units  = units martix (n+2x7)      |
%
% calls: f_co2_atmos
%
% SM, 5/27/2011  
%
function [data,prop,units] = f_add_ttd_Cant(data,prop,units,varargin)

% set default column descritpions ==========================================

sdescr = 'SALNTY';
tdescr = ' THETA';
alk0_descr = '  ALK0';
cant_descr = '  CANT';
co2_sat = 1;

ttdage_descr = 'T12AGE';
ttdwidth_descr = 'T12WID';
ttdageerr_descr = 'AGEERR';
ttdageerl_descr = 'AGEERL';
ttdageerh_descr = 'AGEEHL';
ttdwiderr_descr = 'WIDERR';
ttdwiderl_descr = 'WIDERL';
ttdwiderh_descr = 'WIDEHL';
canterr_descr = 'CANTER';

date_descr = '  DATE';

lat_descr = '   LAT';
lat = [];   % latitude values; can be passed or extracted from 'lat_col' column
hemisphere = 'N';
interp_lat_range = [];

meas_yr = [];     % no measurement year => assumes that column with date exists

include_error = 0;
sol_error_on = 0;
sat_error_on = 0;
atmos_error_on = 0;

%co2_sol_perc_error = ...;
%co2_sat_error =  ...;
%co2_atmos_error = ...;
%ave_flag = 'mean';

age_error_on = 0;
width_error_on = 0;

% check variable input arguments and replace defaults if argument passed ======
i = 1;
while i <= length(varargin)

  if isstr(varargin{i})
    argin = upper(varargin{i});
    if strcmp(argin,'CANT_COL')
      cant_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'ALK0_COL')
      alk0_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'SAL_COL')
      sdescr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'TEMP_COL')
      tdescr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'TTD_MEAN_AGE_COL')
      ttdage_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'TTD_WIDTH_COL')
      ttdwidth_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'SURFACE_SAT')
      co2_sat = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'DATE_COL')
      date_col = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'LAT_COL')
      lat_col = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'LAT')
      lat = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'MEAS_YEAR')
      meas_yr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'INCLUDE_ERROR')
      include_error = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'AGE_ERR_COL')
      ttdageerr_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'AGE_ERL_COL')
      ttdageerl_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'AGE_ERHL_COL')
      ttdageerh_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'WIDTH_ERR_COL')
      ttdwidtherr_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'WIDTH_ERH_COL')
      ttdwidtherh_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'WIDTH_ERL_COL')
      ttdwidtherl_descr = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'SOL_ERROR_ON')
      sol_error_on = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'SOL_PERC_ERROR')
      hin = varargin{i+1};
      co2_sol_perc_error = hin(1);
      sol_error_on = 1;
      i = i+2;
    elseif strcmp(argin,'SAT_ERROR_ON')
      sat_error_on = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'SAT_ERROR')
      hin = varargin{i+1};
      co2_sat_error = hin(1);
      sat_error_on = 1;
      i = i+2;
    elseif strcmp(argin,'ATMOS_ERROR_ON')
      atmos_error_on = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'AVE_UP_DOWN_ERROR')
      ave_flag = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'HEMISPHERE')
      hemisphere = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'INTERP_LAT_RANGE')
      interp_lat_range = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'AVE_UP_DOWN_ERROR')
      ave_flag = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'AGE_ERROR_ON')
      age_error_on = varargin{i+1};
      i = i+2;
    elseif strcmp(argin,'WIDTH_ERROR_ON')
      width_error_on = varargin{i+1};
      i = i+2;
    else
      error_str = ['Unkown identifier ''' varargin{i} '''!'];
      error(error_str);
    end
  else
    error_str = ['Identifier must be a string!'];
    error(error_str);
  end

end


% find columns of ================================================================

[dm,dn] = size(data);
 
ctht=sel_col(prop,tdescr);
csal=sel_col(prop,sdescr);
 
cage=sel_col(prop,ttdage_descr);
cwid=sel_col(prop,ttdwidth_descr);
calk0=sel_col(prop,alk0_descr);

% determine time of mearuserment if no measurement year passed =================
if isempty(meas_yr) 
  cyr = sel_col(prop,date_descr);
  meas_yr = data(:,cyr);
end
if length(meas_yr) == 1
  meas_yr = meas_yr * ones(dm,1);
end

% get atmospheric CO2 concetrations ============================================
% (selection of different hemisphere not on yet)
[pco2_atmos,atmos_yr] = f_co2_atmos;
atmos_n = length(atmos_yr);
[pco2_preind,i] = min(pco2_atmos); % assume that history starts in 1800s
yr_preind = atmos_yr(i);

% interplate too much finer scale ==============================================
% (actually do this later after surface concentrations calculated)
%hatmos_yr = atmos_yr(1):1/100:atmos_yr(end);
%pco2_atmos = sm_interp1(atmos_yr,pco2_atmos,hatmos_yr);
%atmos_yr = hatmos_yr;

% loop over data (follow procedures in f_mk_anthro_gases_ttd_age_lookuptable.m)
clear Cant
for di = 1:dm
% 
  if mod(di,2) == 0
%     di    % just to report some progress

    fprintf('Year %d Point %d/%d aCO2: %f\n',meas_yr(di),di,dm,Cant(di-1))
  end
%   tic;
  % mean age and width for inverse Gaussian
  hgamma = data(di,cage);
  hdelta = data(di,cwid);

  % salinity, temperature, and preformed alkalinity of sample  
  hsal = data(di,csal);
  htemp = data(di,ctht);
  halk0 = data(di,calk0);


  % only do convolution calculations for Cant if all necessary data exists
  if isnan(hgamma) | isnan(hdelta) | isnan(hsal) | isnan(htemp) | isnan(halk0) 

    Cant(di) = NaN;

  else

    % calculate surface time history of anthropogenic CO2
    clear cant_surf
    for i = 1:atmos_n 
      %i
      [hco2,hpco2,hhco3,hco3,dic,halk,hphtotal,hphfree] = f_csys(data(di,ctht),data(di,csal),0,401,pco2_atmos(i)*1e-6,data(di,calk0)*1e-6,0,2); % needs mol/kg; returns mol/kg
      [hco2,hpco2,hhco3,hco3,dic_preind,halk,hphtotal,hphfree] = f_csys(data(di,ctht),data(di,csal),0,401,pco2_preind*1e-6,data(di,calk0)*1e-6,0,2);  % needs molk/kg; returns mol/kg
      cant_surf(i) = (dic - dic_preind)*1e6;  % convert to umol/kg
    end % for i = 1:atmos_n
    
%    f_pause

    % interplate too much finer scale ==============================================
    cant_yr = atmos_yr(1):1/100:atmos_yr(end);
%     cant_surf = sm_interp1(atmos_yr,cant_surf,cant_yr);
    cant_surf = interp1(atmos_yr,cant_surf,cant_yr);
  
    % determine time vector for intergral from 0 to inf. of c(t-tdash), G(tdash) with regard to tdash
    dt = cant_yr(2) - cant_yr(1);  % assume even spacing
    [dummy,ai] = min(abs(cant_yr - meas_yr(di)));  % cut off surface
                        % history at time closest to time of measurement  
                        % (for simplicity don't interpolate; ok since surface 
                        % concentrations interpolated to 1/100 of a year earlier on)
    ai = ai(1);  % in case two values (i.e. meas_yr exactly in middle of two 
                 % surface vals) use first one
    n = ai;  % that's how many dt intervals need to be integrated until 
             % c goes to zero (i.e. don't go all the way to infinity)
    t = cant_yr(ai);
    tdash = [0:1:n]*dt;  % times over which to integrate (n+1 values since 0 
                         % added at end)

    % calculate G and C vectors for convolution integral
    %tic % this step takes on the order of 0.04 seconds (more or less) on laptop
    Gint=sqrt(hgamma^3./(4*pi*hdelta^2*tdash.^3)).*exp(-hgamma*((tdash-hgamma)).^2./(4*(hdelta^2)*tdash));   % G0 in integral; following Waugh et al., (2003), equation 16
    % toc

    Gint(1) = 0; % set first one at tdash=0 to zero since Inverse Gaussian
                 % approaches zero toward t=0, but actually is NaN/unidenti-
                 % fied at t=0 because of t in the denominators (or perhaps
                 % need to do proper limes) 
%     C0int = [sm_interp1(cant_yr,cant_surf,t-tdash)];   % C0 in integral 
    C0int = [interp1(cant_yr,cant_surf,t-tdash)];   % C0 in integral 
                        % (interpolation is a way to pick out values and flip them; 
                        % t-tdash should correspond to times at cant_yr
                        % anyway) 

    if isnan(C0int(end))  % should actually be Nan since it should be outside
                          % the atmos_year since tdash n+1 long;  
      C0int(end) = 0;     % replace last value with zero assuming that tracer
                          % is transient
    end
    C0int = C0int(:);  % make sure that column integral
    Gint = Gint(:);    % make sure that column integral

    % do convolution of integral (integral sum, most simple wat of doing it; same
    % as trapz method
    % tic
    Cant(di) = sum([Gint(1)*C0int(1)*dt/2; 
                            Gint(2:end-1).*C0int(2:end-1)*dt; ...
                            Gint(end)*C0int(end)*dt/2]); % most simple integral
%     toc 
%   fprintf('Year %d Point %d/%d aCO2: %f\n',meas_yr(di),di,dm,Cant(di))
  end  % if isnan(hgamma) | isnan(hdelta) | isnan(hsal) | isnan(htemp) | isnan(halk0) 

end 

% calculate error in Cant  =================================================

dCant = 0*Cant;

% not implemented yet;  see f_mk_anthro_gases_ttd_age_lookuptable.m
% for procedures of error calculations

if include_error
  if sol_error_on

    %df11_sol = f11_sol_perc_error*f11_sol;
    %df12_sol = f12_sol_perc_error*f12_sol;
  %
  %  dpcfc11 = dpcfc11 + (-data(:,cf11)./(f11_sat*f11_sol.^2)).^2 ...
  %                                                       .* df11_sol.^2;
  %  dpcfc12 = dpcfc12 + (-data(:,cf12)./(f12_sat*f12_sol.^2)).^2 ...
  %                                                     .* df12_sol.^2;
  end

  if sat_error_on

  %  dsat11 = f11_sat_error;
  %  dsat12 = f12_sat_error;
  %
  %  dpcfc11 = dpcfc11 + (-data(:,cf11)./(f11_sol.*f11_sat^2)).^2 ...
  %                                                       .* dsat11.^2;
  %  dpcfc12 = dpcfc12 + (-data(:,cf12)./(f12_sol.*f12_sat^2)).^2 ...
  %                                                       .* dsat12.^2;
  end

  if atmos_error_on
  %
  end

end
%dCant = sqrt(dCant);


% add Cant as a columns ======================================================

data(:,dn+1) = Cant; 

if include_error

  %%%  not implemented yet

end
  
% add property and units description of new columns ========================

[mp,np] = size(prop);
pdescr = cant_descr;
nn = length(pdescr);
if nn > np
  prop(mp+1,1:np) = pdescr(1:np); % chop of rest of Cant description
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

if include_error

  pdescr = canterr_descr;
  nn = length(pdescr);
  if nn > np
    prop(mp+2,1:np) = pdescr(1:np); % chop of rest of Cant description
  else
    prop(mp+2,1:nn) = pdescr;
  end

  [mu,nu]=size(units);
  udescr = 'umol/kg';
  nn = length(udescr);
  if nn > nu
    units(mu+2,1:nu) = udescr(1:nu);
  else
    units(mu+2,1:nn) = udescr;
  end

end


