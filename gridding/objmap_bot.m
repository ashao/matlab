function [Grid, Errors, Smooth_Grid, Fine_Grid, opts] = objmap(x, y, z, xgrid, ygrid, opts)
% 10/13/06: worked on implementing error (finally); as indicator to perform 
%           error calculation the option opts.doerror (default 0) is added 
%
% required arguments:
%
% x      -- array of x values (must be organized by station)
% y      -- array of y values
% z      -- array of z values
% xgrid  -- array of xgrid values
% ygrid  -- array of ygrid values
%
% arguments below will default if not passed:
%
% opts.verbose -- 1 = show progress, 0 = don't show
% opts.nseg    -- number of stations in each segment
%
% opts.efold   -- large scale field decay length
% opts.efold2  -- small scale field decay length
% opts.evar    -- large scale field variance (relative to 1)
% opts.evar2   -- small scale field variance (relative to 1)
%
% opts.sdepth  -- standard depths
%
% opts.sdist   -- standard distances (found here from x-array)
% opts.maxdist -- maximum allowed distance to call 1 unit of station spacing
%
% opts.doerror -- 1 = calculate error map, 0=don't do error map
%
%  Roemmich D., Optimal Estimation Of Hydrographic Station Data and Derived Fields,
%  Journal of Physical Oceanography, 13, 1544-1549, Aug 1983.
%

%
% set defaults  (SM, 5/23/07: leave them as is although efold2=1, efold=2 and
% different evar (originally evar2 =0.3) might make more sense -> defaults need
% to be changed from driver routine)
%
out = 0;
nseg = 12;
efold = 40;  % correlation scales should be 40*standard depth/dist interval
efold2 = 2;  % correlation scales should be 2*standard depth/dist interval
%efold2 = 1;   % make all information on correlation scales be in sdist and sdepth (i.e. division by 1 when calculating covariance does not change anything)
evar = 0.1;
evar2 = 0.1;
epower = 1;
epower2 = 1;
doerror = 0;
%
maxdist = 200;	% this assumed distance and km  (versus lat, lon)
%
sdepth = [0 50 100 150 200 250 300 350 400 450 500 550 600 700 800 900 1000 1100  1200 1400 1600 1800 2000 2200 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000];

if (nargin < 6)
   opts = [];
end

if isfield(opts, 'verbose'), out = opts.verbose; end
if isfield(opts, 'efold'),  efold  = opts.efold; end
if isfield(opts, 'efold2'), efold2 = opts.efold2; end
if isfield(opts, 'evar'),   evar   = opts.evar; end
if isfield(opts, 'evar2'),  evar2  = opts.evar2; end
if isfield(opts, 'epower'), epower = opts.epower; end
if isfield(opts, 'epower2'), epower2 = opts.epower2; end
if isfield(opts, 'maxdist'), maxdist = opts.maxdist; end
if isfield(opts, 'sdepth'),  sdepth = opts.sdepth; end
% added:
if isfield(opts, 'nseg'),  nseg = opts.nseg; end
if isfield(opts, 'doerror'),  doerror = opts.doerror; end

efold
efold2
evar
evar2
epower
epower2
maxdist
sdepth
nseg
doerror

if (~out)
   out = fopen('objmap.out', 'w');
end

%
% ensure all are row vectors
%
x = x(:)'; y = y(:)'; z = z(:)';
xgrid = xgrid(:)'; ygrid = ygrid(:)';
%

n = length(x);
ngrid = length(xgrid);

% 8/3/99: possible to pass field sdist now
if isfield(opts, 'sdist')
   sdist = opts.sdist;
   idist = [1:length(sdist)];
   sta_index = find(abs(diff(x)) > 1e-5) + 1;
   sta_index = [1 sta_index];
   nsta = length(sta_index);
else
  %
  % find distance between stations (sdist array)
  % (should be passed, since for CFCs not all stations have samples)
  %
  nsta = 1;
  sdist(1) = x(1);
  idist(1) = 1;
  lastx = x(1);
  sta_index(1) = 1;
  for i = 2: n
    if (x(i) ~= lastx)
       nsta = nsta + 1;
       sdist(nsta) = x(i);
       sta_index(nsta) = i;
       %
       % adjust standard station spacing index
       %
       if (x(i) - lastx < maxdist)
          idist(nsta) = idist(nsta-1) + 1;
       else
          idist(nsta) = idist(nsta-1) + (x(i) - lastx)/maxdist;
       end

       lastx = x(i);
    end
  end
end

sdist
idist
sta_index
length(sta_index)
nsta

sta_index(nsta + 1) = n + 1;

%
fprintf(1, 'Gridding %d stations with %d points\n', nsta, n);
%
%
% Convert x, y, xgrid and ygrid to units of standard spacing
%
fprintf(out, 'Converting x and y into units of station spacing and standard depth\n');
%
x = interp1(sdist, idist, x);
y = interp1(sdepth, [1:length(sdepth)], y);
%
fprintf(out, 'Converting xgrid and ygrid into units of station spacing and standard depth\n');
%
xgrid = interp1(sdist, idist, xgrid);
ygrid = interp1(sdepth, [1:length(sdepth)], ygrid);
%
% create blank grids
%
Grid = zeros(length(ygrid),length(xgrid));
if (nargout > 1)
   Errors = zeros(length(ygrid),length(xgrid));
end
if (nargout > 2)
   % original:
   %Smooth_grid = zeros(length(ygrid),length(xgrid));
   % changed
   Smooth_Grid = zeros(length(ygrid),length(xgrid));
end
% added:
if (nargout > 3)
   Fine_Grid = zeros(length(ygrid),length(xgrid));
end

% 
% Break stations into segments here
%
Nsegment = nsta - nseg + 1;
if (Nsegment <= 0)
   Nsegment = 1;
end
for sta_start=1: Nsegment

    sta_end = sta_start + nseg - 1;
    if (sta_end > nsta)
       sta_end = nsta;
    end
    istart = sta_index(sta_start);
    iend = sta_index(sta_end + 1) - 1;
%
% capital letters are for this segment
%
    X = x(istart : iend);
    Y = y(istart : iend);
    Z = z(istart : iend);
    N = iend - istart + 1;
%
% note x is now in units of station index (units of station spacing)
%
    old_version = 0;
    if old_version
    grid_start = sta_start + nseg/2 - 1;
    grid_end = grid_start + 1;
    if (sta_start == 1)
       grid_start = 1;
    end

    if (sta_end == nsta)
       grid_end = nsta;
       Xgrid = find(xgrid >= grid_start & xgrid <= grid_end);
    else
       Xgrid = find(xgrid >= grid_start & xgrid < grid_end);
    end
    end

    % determine x-grid within center of data window
    %
    grid_start = [];
    grid_end   = [];

    if (sta_start == 1)
      if min(xgrid) < x(sta_index(1))
        grid_start = min(xgrid); % include grid points outside station domain
      else
       grid_start = x(sta_index(1));
      end
    end

    if (sta_end == nsta)
      if max(xgrid) > x(sta_index(nsta))
        grid_end = max(xgrid);  % include grid points outside station domain
      else
        grid_end = x(sta_index(nsta));
      end
    end

    if 2*floor(nseg/2)==nseg     % even numbers of segments
      wi1 = sta_start + nseg/2 - 1;
      if isempty(grid_start)
        grid_start = x(sta_index(wi1));
      end
      if isempty(grid_end)
        grid_end = x(sta_index(wi1+1));
      end
    else                           % uneven numbers of segments
      wmi  = sta_start + nseg/2 - 0.5;
      if isempty(grid_start)
        grid_start = mean(x(sta_index(wmi -1 : wmi)));
      end
      if isempty(grid_end)
        grid_end = mean(x(sta_index(wmi : wmi+1)));
      end
    end

    if (sta_end == nsta)
      Xgrid = find(xgrid >= grid_start & xgrid <= grid_end);
    else
      Xgrid = find(xgrid >= grid_start & xgrid < grid_end);
    end

    if (~isempty(Xgrid))
%
       fprintf(out, 'Using data from stations %d and %d (%d points)\n', sta_start, sta_end, N);
       fprintf(out, 'To evaluage grid between stations %d and %d\n', grid_start, grid_end);
%
       Acov = zeros(N,N);
       for i = 1: N
           dist = sqrt((X(i) - X).^2 + (Y(i) - Y).^2);
           Acov(i,:) = exp(-(dist/efold).^epower);
           Acov(i,i) = Acov(i,i) + evar;
       end
%
%
% calculate spatially weighted mean and subtract it from the data
%
       fprintf(out, 'Calculating spatially weighted mean: ');
       W = Acov\Z';
       Wsum = sum(W);
       Asum = sum(Acov\ones(N, 1));
       Mean = Wsum/Asum;
%
       fprintf(out, 'Wsum = %f, Asum = %f, Mean = %f\n', Wsum, Asum, Mean);
%
       Z = Z - Mean;
%
%
% get large scale field weights
%
       fprintf(out, 'Determining large-scale field weights\n');
       W = (Acov\Z');
%
% subtract smooth field from data
%
       C = zeros(1,N);
       for i = 1:N
           dist = sqrt((X(i) - X).^2 + (Y(i) - Y).^2);
           C = exp(-(dist/efold).^epower);
           Z(i) = Z(i) - C*W;
       end
%
% get small scale field weights
%
       for i = 1: N
           dist = sqrt((X(i) - X).^2 + (Y(i) - Y).^2);
           Acov(i,:) = exp(-(dist/efold2).^epower2);
           Acov(i,i) = Acov(i,i) + evar2;
       end
%
       fprintf(out, 'Determining small-scale field weights\n');
       % W2 = (Acov\Z');
       % calculate explicit inverse instead since also used for error
       Ainv = Acov\eye(size(Acov));
       W2 = Ainv*Z';
%
       if doerror    
         % caclulate some variables used for estimating error for each
         % grid point (which according to Roemmich 1983, p. 1546 is "swept 
         % into the small-scale component of the field"), SM 10/13/06
         % Ainv = Acov\eye(size(Acov)); %do above
         Ainv2 = Ainv*ones(N,1);  % same as Acov\ones(N,1);
         Asum = sum(Ainv2);
       end
       fprintf(out, 'Evaluating xgrid at indices: ');
       fprintf(out, '%d ', Xgrid);
%
% evaluate at xgrid, ygrid
%
       ydist = zeros(length(ygrid), length(Y));
       for j = 1: length(ygrid)
           ydist(j,:) = (ygrid(j) - Y).^2;
       end
%
       C2 = zeros(1,N);
       for i = Xgrid
           xdist = (xgrid(i) - X).^2;
%
           for j = 1: length(ygrid)
               dist = sqrt(xdist + ydist(j,:));
%
               if (epower == 1)
                  C = exp(-(dist/efold));
               else
                  C = exp(-(dist/efold).^epower);
               end
%
               if (epower2 == 1)
                  C2 = exp(-(dist/efold2));
               else
                  C2 = exp(-(dist/efold2).^epower2);
               end
%
% 8/9/99: changed to sum of smooth and fine grid (cf. below) because of spped
     %          Grid(j,i) = Mean + C*W + C2*W2;
% original
% Could return Grid (smooth) and Grid2 (fine) grids
%
%              Grid(j,i) = Mean + C*W;
%              Grid2(j,i) = C2*W2;
% changed:
% DO return Grid (smooth) and Grid2 (fine) grids:
%
               Smooth_Grid(j,i) = Mean + C*W;
               Fine_Grid(j,i) = C2*W2;
               Grid(j,i) = Smooth_Grid(j,i) + Fine_Grid(j,i);

               if doerror
                 % following (2) in Roemmich (1983)
                 %Errors(j,i) = 1-C2*(Acov\C2')+(1-C2*Ainv2)^2/Asum; 
                 Errors(j,i) = 1-C2*(Ainv*C2')+(1-C2*Ainv2)^2/Asum; 
                 % multiply by standard deviation to get error in units
                 % acc. to John Lyman (see e.g. Wunsch 1996 book, p. 116,
                 % equ. 3.3.11); use e.g. spatial stdv on depth surfaces (!?);
                 % however percent error (of variance) which it is (as is) is 
                 % in some ways more meaningful; when adding errors could
                 % add them in quadrature and then devide by stdv's/variance
               end
           end
       end

       fprintf(2, 'Completed %4.1f percent\n', 100*Xgrid(length(Xgrid))/ngrid);

   end	% isempty

end % segment

if nargout > 4  % return some of the options    (added 05/08/98)
  opts.efold = efold;
  opts.efold2 = efold2; 
  opts.evar = evar; 
  opts.evar2 = evar2; 
  opts.epower = epower;
  opts.epower2 = epower2;
  opts.maxdist = maxdist;
  opts.sdepth = sdepth;
  opts.sdist = sdist;
end
