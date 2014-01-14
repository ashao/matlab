%
% JLAB_CHANGES   Changes to JLAB in each release.
%   _________________________________________________________________
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
%   Changes new in version 0.93 (current release)
%
%   Version 0.93 includes support for three new published or in press papers, 
%   as well as some minor changes.
%
%   New functions
%
%   makefigs_vortex        - Makes figures for Lilly, Scott, and Olhede (2011), GRL.
%   makefigs_multivariate  - Makes figures for Lilly and Olhede (2012), ITSP.
%   makefigs_trivariate    - Makes figures for Lilly (2011), ITSP.
%   ridgelen      - Wavelet ridge length expressed as number of full cycles.
%   periodindex   - Returns time index in increments of instantaneous period.
%
%   Minor changes and improvements
%
%   READ_SAND updated for version 14.1 of Smith and Sandwell dataset.
%   CELLPLOT rewritten for speed improvement for large datasets.
%   INSTFREQ fixed to improve performance under periodic boundary conditions.
%   CELLPLOT fixed bug for longitude shifting option.   
%   JLAB_RUNTESTS fixed bug causing test-running script to fail on some systems.
%   JMAT3 fixed inconsistent sign convection for rotation about y-axis.
%   CELL2COL now enforces convention of using INFs (not NANs) for missing data.
%   JOINTFREQ added test, fixed incorrect documentation for bandwidth.
%   MSPEC increased number of iterations to address some convergence issues. 
%
%  _______________________________________________________________________
%
%   Changes new in version 0.92 
%
%   Version 0.92 is a major new release, including:
% 
%       -- Complete online help files at http://www.jmlilly.net/doc/.  
%       -- A ``highlights'' file, JLAB_HIGHLIGHTS, introducing the best of JLAB
%       -- Enhancements to wavelet and signal analysis software
%       -- New tools for treating bivariate and trivariate time series in JELLIPSE
%       -- Mapping scattered data in 2D or on the sphere with POLYSMOOTH
%       -- New module, JCELL, for cell arrays of numeric data
%       -- New functions and significant bugfixes
%       -- Scripts for making figures from four new papers:
%
%   makefigs_multivariate  - Makes figures for Lilly and Olhede (2012), ITSP.
%   makefigs_trivariate    - Makes figures for Lilly (2011), ITSP.
%   makefigs_analytic      - Makes figures for Lilly and Olhede (2010b), ITIT.
%   makefigs_bandwidth     - Makes figures for Lilly and Olhede (2010a), ITSP.
%
%   New functions
%
%   jlab_highlights - Meet the best rountines in JLAB.  You really want to read this.
%   lininterp  - Fast linear interpolation for arbitrary-sized arrays.
%   jointfreq  - Joint instantaneous frequency, bandwidth, and curvature.
%   ellsig     - Creates a modulated elliptical signal in two or three dimensions.
%   ellband    - Bandwidth of modulated elliptical signals in two or three dimensions.
%   ellparams  - Ellipse parameters of a modulated bivariate or trivariate oscillation.
%   normvect   - Unit normal vector to the ellipse plane in three dimensions.
%   growl      - Display text message with the Growl notification system.
%   fourier    - The Fourier frequencies for a given length time series.
%   twodstats  - Mean, variance, and covariance of functions of two variables.
%
%   Minor changes and improvements
%
%   TWODMEAN functionality extended to covariance analysis, and renamed TWODSTATS.
%   XYZ2LATLON modified to work with vector of any length.
%   DENSCONT now plots the freezing line of seawater.
%   SPHEREDIST now uses Haversine formula for improved precision at small distances.
%   ANATRANS fixed for nonzero-mean time series input.
%   LATLON2XY and LATLON2UV now support complex-valued input and output.
%   LATLON2UV improvements and added tests.  Refined central difference.     
%   SLIDETRANS rewritten, added extra functionality.
%   SLEPTAP now supports bandpass normalization.
%   MSPEC rewritten, added sample rate normalization, and added new tests.
%   CELLPLOT fixed bug for cell array of length one, compatibility with M_MAP.
%   WAVETRANS made initialization of missing data to complex NANs consistent.
%   TWODHIST, TWODMEAN, TWODMED, added tests for data points outside bins edges.
%   TWODMEAN major bugfix for data with points outside of bin edges.
%   TWODMED bugfix for input row vectors.
%   ELLIPSEPLOT increased plot speed for input style arguments.
%   INSTFREQ now handles joint moments of multivariate signals.
%   RIDGEMAP initialization consistency fixed.
%   LINESTYLE bugfix for some input symbol types
%   RIDGEMAP change to output format for joint ridges.
%   DAWSON fixed non-executing tests.
%   ASCII2NUM bugfix for multiple input arguments.
%   CATSTRUCT extended to multiple input arguments.
%   SPHEREDIST now also works for computing along-track distance.
%   VZEROS now has option to initialize to INFs.
%   XOFFSET, YOFFSET line offset order reversed so first line doesn't move.
%   ELLBAND bugfix for missing DT in time derivatives, and new tests.
%   ELLIPSEPLOT now automatically sets aspect ratio.
%   ANATRANS added tests and flag for mirror and other boundary conditions.
%   ELLSIG added test for inversion of ELLPARAMS in trivariate case. 
%   ISCOL and ISROW renamed JISCOL and JISROW.
%  _______________________________________________________________________
%
%   Changes new in version 0.91 
%
%   Version 0.91 is a major new release, including:
%
%       -- Significantly enhanced wavelet and ridge analysis software
%       -- New tools for treating bivariate time series in JELLIPSE
%       -- Mapping scattered data in 2D or on the sphere with POLYSMOOTH
%       -- New module, JCELL, for cell arrays of numeric data
%       -- Dozens of new functions and significant fixes
%
%   *Important Notes*
%
%   --- All wavelet functions that previously used cyclic frequency, 
%   as in cos (2 pi f t), now use radian frequency, as in cos (omega t).
%
%   --- ELLDIFF, ELLRAD, ELLVEL, and ELLIPSEPLOT, and ELLSIG now all use 
%   KAPPA and LAMBDA to specify the ellipses rather than A and B.
%
%   Scripts for making figures from new papers
%
%   makefigs_analytic   - Makes figures for Lilly and Olhede (2010c), ITIT.
%   makefigs_asilomar   - Makes figures for Lilly and Olhede (2010b), ACSS.
%   makefigs_bandwidth  - Makes figures for Lilly and Olhede (2010a), ITSP.
%   makefigs_morsies    - Makes figures for Lilly and Olhede (2009), ITSP.
%
%   New functions -- from JARRAY
%
%   matinv     - Fast inversion of arrays of small sub-matrices.
%   matmult    - Matrix multiplication for arrays of sub-matrices.
%   vmedian    - Median over finite elements along a specified dimension.
%
%   New functions -- from JSTATS
%
%   polysmooth - Smoothing scattered 2D data with local polynomial fitting.
%   bindata    - Rapidly sort data into adjacent bins.
%   twodsort   - Distances from data points to nearby grid points.
%
%   New functions -- from JPLOT
%
%   ztick      - Set locations of z-axis tick marks.
%   dlines     - Add diagonal lines to a plot.
%
%   New functions -- from JSPHERE
%
%   spheresort - Sorted great circle distances to nearby points on a sphere.
%   spherelap  - Laplacian of a field on the surface of a sphere.
%   spherecurl - Curl of a vector field on the surface of a sphere.
%   cellplot   - Plot elements of a cell array of numeric arrays.
%
%   New functions -- from JOCEANS
%
%   corfreq    - Coriolis frequency in cycles per hour.
%   ssh2eke    - Converts alongtrack sea surface height to eddy kinetic energy.
%   read_sand  - Read topography data from the Smith and Sandwell Database.
%
%   New functions -- from JSIGNAL
%
%   morlfreq   - Compute Morlet wavelet carrier frequency given peak frequency.
%   morsehigh  - High-frequency cutoff of the generalized Morse wavelets.
%   powermean  - Power-weighted along a specified dimension.
%   morsespace - Logarithmically-spaced frequencies for generalized Morse wavelets.
%
%   New module --- JCELL, for operating on cell arrays of numeric arrays.
%
%   col2cell   - Converts 'column-appended' data into cell arrays of column vectors.
%   cell2col   - Converts cell arrays of column vectors into 'column-appended' data.
%   cellaresame - Tests whether two cell arrays of arrays are the same.
%   cellength  - Length of each element in a cell array.
%   cellsize   - Size of each element in a cell array along specified dimension.
%   cellmax    - Maximum of each element in a cell array.
%   cellmin    - Minimum of each element in a cell array.
%   cellreal   - Real part of each element in a cell array.
%   cellimag   - Imaginary part of each element in a cell array.
%   cellabs    - Absolute value of each element in a cell array.
%   cellangle  - Complex argument (angle) of each element in a cell array.
%   cellconj   - Complex conjugate of each element in a cell array.
%   cellplot   - Plot elements of a cell array of numeric arrays.
%   
%   Minor changes and improvements
% 
%   READ_PATHFINDER region definition change for consistency.
%   ELLRAD and ELLVEL fix for cell arrays.
%   JLAB_RUNTESTS improvements.
%   DEG180, DEG360 speed improvements.
%   WAVETRANS fixed bug that rotated transform with real-valued wavelets by sqrt(-1).
%   PDFPROPS major speed increase.
%   TWODMEAN and TWODMED bugfix for data outside of bins.
%   TWODMEAN bugfix for fast algorithm on account of annoying CUMSUM misbehavior.
%   SPHERECURL and SPHEREDIV now output components, e.g. for computing strain.
%   Replaced incorrect usage of FINDSTR with STRFIND throughout. 
%   COL2MAT and MAT2COL now return empty given empty.
%   MORSEMOM now also computes energy cumulants.
%   INSTFREQ now outputs higher-order instantaneous modulation functions.
%   WAVETRANS added change in normalization for complex-valued time series.
%   MATMULT speed improvement.
%   INREGION region vector definition changed to match Matlab's AXIS.
%   MORSEWAVE and LAGUERRE improved to handle very small values of GAMMA.
%   TWODMEAN, TWODMED, and TWODHIST speed increase by stripping non-finite data.
%   TWODMEAN, TWODMED, and TWODHIST now have flag for memory-efficient algorithm.
%   VMEDIAN bugfix for columns or rows of all NaNs.
%   NUM2YF renamed YEARFRAC.
%   FINDFILES now has option for include or excluding specific strings.
%   ELLIPSEPLOT style handling and argument sorting improvements.
%   VECTMULT now handles all dimensions, and speed improvements.
%   ELLRAD and ELLVEL output argument order change.
%   VCELLCAT now works with some empty sets.
%   LINEHANDLES, LINESTYLE, and LINERING bugfixes.
%   LINEHANDLES, LINESTYLE, and LINERING now work with contours (patch objects).
%   SPHEREDIV and SPHEREGRAD now have options for handling boundary points.
%   RAD2DEG renamed JRAD2DEG to avoid conflict with Matlab's mapping toolbox.
%   DEG2RAD renamed JDEG2RAD to avoid conflict with Matlab's mapping toolbox.
%   BELLBAND renamed MODFUN.
%   JINTERP fixed bug for matrix input; this fixed also PDFMULT.
%   VINDEX now returns empty output for empty input.
%   ELLVEL input argument change for optional leading DT.
%   VDIFF input argument change for optional leading DT.
%   MORSEWAVE change in phase convention for wavelets higher than zeroth order.
%   NOXLABELS and NOYLABELS improved to allow empty max or min label.
%   LENGTHCELLS renamed CELLENGTH.
%   VDIFF modified for variable behavior at endpoints.
%   LATLON2UV modified to prevent appearance of NANs at first and last times.
%   SLEPTAP refactoring, tests, and support for cell output.
%   SIMPLEPDF now supports exponential PDF.
%   TWODMEAN refactoring and additional speed and memory usage improvements.
%   MORLWAVE improvements; note change to input arguments.
%   MORLWAVE support for negative frequencies, and tests, added.
%   MORSEWAVE support for negative frequencies, and tests, added.
%   MORSEWAVE support for non-unit sample rate removed.
%   MORSEBOX output arguments set to NAN when imaginary.
%   PACKROWS, PACKCOLS, PACKBOTH: fixed bug switching box back on.
%   MATMULT substantially improved and changed:
%      Now support multiplication of non-square matrices of any compatible size. 
%   RIDGEWALK bugfix for occasional sign reversals in frequency.
%   RIDGEWALK now no longer expecting chaining parameter ``ALPHA'' as input.
%   RIDGEWALK refactored and simplified.
%   RIDGEWALK behavior for 3-D transforms has now changed.
%   RIDGEWALK now computes joint ridges.
%   RIDGEINTERP no longer interpolates transform W by default. 
%   VNAN added check for all input arrays being the same size.
%   ELLIPSEPLOT change to aspect ratio argument for more flexibility.
%   ELLIPSEPLOT order of magnitude speed increase for multi-column data.
%   RIDGEMAP simplication and input argument change.
%   LINESTYLE bug for line handle input fixed.
%   WAVESPECPLOT contouring improved to fix Matlab contour closing bug.
%   INSTFREQ now outputs frequency, bandwidth, and chirp rate.
%   INSTFREQ now no longer supports anti-analytic signals.
%   RIDGEWALK bugfix in test case.
%   SIMPLEPDF bugfixes for non-standard Cauchy distribution.
%   ELLDIFF simplifications.
%
%   Removed functions
%
%   INDEXAND, INDEXOR, ISBOOLEAN removed.
%   MORSEREGION removed.
%   SLEPWAVE removed; use MORSEWAVE instead.
%   UNWRANGLE removed; use Matlab's UNWRAP.
%   FINDFIRST and FINDLAST removed; use Matlab's FIND with 'first' or 'last'.
%   ELLCONV removed; use simpler ELLPARAMS.
%   WALPHA removed.
%   WAVERECON removed.
%   NORMFORM, TRANSCONV, MANDN, RANDSPECMAT removed.
%   VECTMULT3 removed; use VECTMULT.
%   FONTSIZESETS removed; incorportated into JLAB_SETTINGS.
%   MODMAX, MODMAXPEAKS, ISOMAX, and EDGEPOINTS removed.
%   ISRIDGEPOINT removed; absorbed into RIDGEWALK.
%   RIDGEDEBIAS removed; absorbed into RIDGEWALK.
%   RIDGESTRUCT removed; absorbed into RIDGEWALK.
%   PATCHHANDLES removed; now implemented by LINEHANDLES.
%   SECONDAXIS removed.
%   NDREP removed; use VREP.
%   ELLRIDGE removed; this functionality is now carried out by RIDGEWALK.
%   NEARESTPOINT removed.
%   BANDNORM removed; absorbed into MORSEWAVE.
%   MODFUN removed; this functionality is now carried out by INSTFREQ.
%   _______________________________________________________________________
%
%   Changes new in version 0.90 (Previous release)
%
%   Version 0.90 is a major new release, including:
%
%       -- Improved organization of modules
%       -- Dozens of new functions and significant fixes
%       -- New functions of spherical geometry
%       -- Fast interpolation functions: QUADINTERP and CUBEINTERP
%       -- Refactoring and improvements to wavelet ridge code
%       -- Additional Morse wavelet functions
%
%   New functions --- Generalized Morse Wavelets
%                        In support of Lilly and Olhede (2009)
%
%   morsemom      - Frequency-domain moments of generalized Morse wavelets.
%   morsederiv    - Frequency-domain derivatives of generalized Morse wavelets.
%   morsexpand    - Generalized Morse wavelets via a Taylor-series expansion.
%   morsebox      - Heisenberg time-frequency box for generalized Morse wavelets.
%   morseprops    - Properties of the demodulated generalized Morse wavelets.
%   morsefreq     - Frequency measures for generalized Morse wavelets.  
%   dawson        - The Dawson function. [By P. J. Acklam]
%   dawsonderiv   - Derivatives of the Dawson function.
%   makefigs_morsies - Makes figures for Lilly and Olhede (2009).
%
%   New functions  -- Instantaneous frequency and bandwidth
%                        In support of Lilly and Olhede (2010a)
%
%   instfreq    - Instantaneous frequency and bandwidth of an analytic signal.
%   bellpoly    - Complete Bell polynomials.
%   makefigs_analytic - Makes figures for Lilly and Olhede (2010a).
%
%   New functions --- Spherical geometry
%
%   xyz2latlon - Converts Cartesian coordinates into latitude and longitude.
%   latlon2xyz - Converts latitude and longitude into Cartesian coordinates.
%   uvw2sphere - Converts a 3D Cartesian vector to a 3D spherical vector.
%   sphere2uvw - Converts a 3D spherical vector to a 3D Cartesian vector.
%   uvw2hor    - Projects a 3D Cartesian vector into a horizontal vector on a sphere.
%   hor2uvw    - Converts a horizontal vector on a sphere into a 3D Cartesian vector.
%   spherediv  - Divergence of a vector field on the surface of a sphere.
%   spheregrad - Gradient of a field on the surface of a sphere.
%   neareastpoint - Finds the nearest point to a specified point on the sphere.
%
%   New functions --- Aquarius satellite
%
%   aquaplot   - Plot Aquarius satellite radiometer footprint.
%   aquaprint  - Compute Aquarius satellite radiometer footprints.
%   aquasal    - Aquarius salinity change with brightness temperature.
%
%   New functions
%
%   lonshift   - Shifts longitude origin for plotting purposes.
%   heat2evap  - Transform latent heat loss into units of evaporation.
%   ellridge   - Extract "elliptical" ridges for bivariate time series.
%   findfiles  - Returns all files in a directory with a specified extension.
%   quadinterp - Fast quadratic interpolation for arbitrary-sized arrays.
%   cubeinterp - Fast cubic interpolation for arbitrary-sized arrays.
%   ellband    - Elliptical bandwidth of bivariate signal or wavelet transform.
%   jmat3      - 3-D rotation matrix through specified angle.
%   vectmult3  - Matrix multiplication for arrays of three-vectors.
%   ab2kl      - Converts A and B to ellipse parameters Kappa and Lambda.
%   kl2ab      - Converts ellipse parameters Kappa and Lambda to A and B.
%   twodmed    - Median value of a function of two variables.
%   ellsig     - Creates an elliptical signal from ellipse parameters.
%   choose     - Binomial coefficient: CHOOSE(N,K) = "N choose K"
%   latratio   - Set plot aspect ratio for latitude / longitude plot.
%   mom2cum    - Convert moments to cumulants.
%   cum2mom    - Convert cumulants to moments.
%
%   New low-level signal analysis functions
%
%   ridgequantity - Returns the quantity to be minimized for ridge analysis.
%   randspecmat   - Generates random 2x2 spectral matrices for testing
%
%   Minor changes and improvements
%
%   DISCRETECOLORBAR bugfix for Matlab 7.5.
%   JLAB_RUNTESTS now reports summary statistics
%   MORSEWAVE now uses 'bandpass' normalization by default
%   MORSEWAVE support of non-unit sample rate and zero beta
%   MORSEWAVE additional testing with analytic expressions
%   MORSEWAVE normalization flag added
%   HERMPOLY, HERMFUN, HERMEIG changed to output N+1 terms 
%   WAVETRANS now has direct computation of generalized Morse wavelets
%   JDEG2RAD and JRAD2DEG modified to preserve NANs and INFs
%   MORSEFREQ output argument change; new frequency measures now computed   
%   WIGDIST modified to allow for negative frequencies; new sample figure
%   RIDGEWALK now has improved ridge control parameters
%   RIDGEWALK modification --- 
%         Connection of ambiguous ridges now done by minimizing
%         log-transform curvature, rather than amplitude difference
%   RIDGEWALK modification --- RIDGEINTERP now called internally by default
%   RIDGEINTERP now uses superior fast quadratic interpolation via QUADINTERP
%   RIDESTRUCT and RIDGEWALK use improved defintion of ridge "length"
%   RIDESTRUCT and RIDGEWALK use improved chaining algorithm 
%   RIDGEWALK and RIDGEINTERP complete retooling into low-level components
%   RIDGEWALK output argument change ---
%         Transform frequency no longer output; use INSTFREQ. 
%   RIDGEWALK recommenting, added test
%   RIDGEINTERP reclassified as "low-level" function
%   RIDGEPRUNE removed --- use RIDGESTRUCT for same functionality
%   RIDGEWALK bugfix for chaining very sparse ridges (due to VINDEX bug)
%   RIDGEWALK ridge structure format changed for smaller size
%   RIDGEMAP now supports elliptical ridge structure from ELLRIDGE
%   VINDEX bugfix for empty indicies
%   ELLIPSEPLOT now draws ellipses in default colors
%   ELLIPSEPLOT now decimates when input a row index array
%   PACKCOLS bugfix
%   PDFMULT bugfix
%   ELLCONV factor of 2 bugfix for two-parameter rotary form
%   ELLCONV test improvements
%   ELLIPSEPLOT added 'npoints' option to set number of points in ellipse
%   ELLIPSEPLOT improved handling of string input arguments
%   ELLIPSEPLOT supress plotting of missing or zero-radius ellipses
%   BLOCKLEN added block number NUM output and tests 
%   VECTMULT input argument convention change and improved flexibility
%   VFILT added options for mirror and periodic boundary conditions
%   DEG180 and DEG360 modified to accept any input degree range
%
%   Depricated functions
%
%   SLEPENVWAVE
%   TRACK2GRID
%   RIDGERECON
%   _______________________________________________________________________
%
%   Changes new in version 0.85 
%
%   New functions
%
%   fastcontour - Lightning-fast "fake" contouring for large matrices.
%   whichdir    - Returns directory name containing file in search path.
%   range       - RANGE(x)=[MIN(x(:)) - MAX(X(:))];
%   inregion    - Tests whether lat/lon points lie within a specified box.
%   turningpoint  - True for turning points, i.e. local extrema, along rows.
%   crossings     - Find values of an array crossing a specified threshold.
%   orbitbreaks   - Separate orbit into passes based on turning points.
%   to_grab_from_caller - Returns a string to grab variable values from caller.
%   deg360     - Converts degrees from the range [-180,180] to [0,360].
%   deg180     - Converts degrees from the range [0,360] to [-180,180].
%   jrad2deg    - Converts radians to degrees.
%   jdeg2rad    - Converts degrees to radians.
%   latlon2xyz - Convert latitude and longitude into Cartesian coordinates.
%   xyz2latlon - Convert Cartesian coordinates into latitude and longitude.
%   latlon2zeaz - Compute zenith and azimuth angles for satellite beam.
%   zeaz2latlon - Compute latitude and longitude viewed by satellite beam.
%   ascii2num  - Convert ASCII values for numbers into numeric values.
%   jlab_settings - Specifies settings for customizable JLAB properties.
%
%   Major bugs from version 0.84
%
%   LATLON2UV 0.84 was incorrect; corrected & tested in version 0.85.
%   TRACKFILL 0.84 was incorrect; temporarily removed.
%
%   Minor changes and improvements
%
%   Global refactoring based on MLINT suggestions
%   'jlab_runtests figures' now makes all possible sample figures
%   VINTERP renamed to JINTERP to avoid filename clash
%   VINDEX changed to suppress indexing along singleton dimensions
%   STRS2SRAY rewritten to correctly format text containing commas
%   Packing turned off in STRS2SRAY, STRS2MAT, STRS2LIST
%   MAKE modified to accept cell array input
%   LAT/LON function improvements:
%      XY2LATLON and LATLON2XY now use full spherical geometry
%      LATLON2UV rewritten to account for full sphereical geometry
%   COMMENTLINES changed to infer m-files
%   XY2ELAZ and ELAZ2XY replaced with LATLON2ZEAZ and ZEAZ2LATLON
%   EL2DIST and EL2INC renamed ZE2DIST and ZE2INC
%   FILLBAD bugfix for one bad data block; added tests
%   VFILT changed to support N-D matrices
%   VFILT input argument change -- 
%        Use VFILT(X,F,'NONANS') instead of VFILT(X,F,1) 
%   VINT2STR renamed VNUM2STR
%   VSWAP now treats +INF and -INF separately
%   JLAB_RUNTESTS refactored
%   AG2BC test added  
%   LINESTYLE now supports handle input with '-h' flag
%   LINESTYLESETS and COLORSETS now incorporated into JSETTINGS
%   JCONTOUR, JCONTOURF, LINEHANDLES, PATCHHANDLES, all modified
%       to account for changes new in Matlab 7
%   MAKEFIGS_LABCONV bugfixes
%   
%   Depricated functions
%
%   JARROW (available as "ARROW" by R.S. Oldaker, online) 
%   ISSCAL removed; provided by built-in ISSCALAR
%   SPHEREPROJ removed; functionality wrapped into LATLON2XY
%   MAKEFIGS_OLSEF removed pending future improvements
%   TRACKFILL removed pending future improvements
%   _______________________________________________________________________
%
%   Changes new in version 0.84 
%
%   MAKEFIGS_RIDGES released; creates figures for Lilly and Gascard 2006
%   New functions for wave triad interactions in JOCEANS:
%       OM, KMIN, GC_PARAMS, TRIADRES, VTRIADRES, ISRES, RESCOEFF,
%       TRIADEVOLVE, KUN, HFUN, DFUN, I2SS, SS2I, DMSPEC, DMSTD, 
%       DMSKEW, DMASYM
%   WAVETRANS added output size test
%   WAVETRANS output change -- K and M dimensions swapped
%   MSVD added output size tests
%   MSVD input and output change -- K and M dimensions swapped
%   MSVD refactoring
%   MSVD added 'quiet' option
%   ANATRANS now different for real and complex input, as in LG06
%   TESTSERIES added replicated-modulated chirp signal
%   BANDNORM bugfix for high-frequency spillover
%   _______________________________________________________________________
%
%   Changes new in version 0.83
%
%   New function: RIDGEMAP
%   New function: ISLARGEST
%   New functions: IMLOG, RELOG, UNWRANGLE
%   New function: RIDGEPRUNE
%   MAKELLIPSE removed
%   ECCONV bugfix for zero eccentricity signals
%   RIDGEINTERP output argument changes
%   CIRC2ELL absorbed into ELLIPSEPLOT
%   PF_PARAMS changed to output |LON|<180
%   VDIFF timestep functionality added
%   RIDGEWALK, RIDGEINTERP, RIDGEMAP modified to work with ridge structures 
%   ISMAT bugfix
%   TWODHIST bugfix for negative data
%   HILTRANS changed to support matrix input
%   VSUM bugfix for no NaNs 
%   STRCAPTURE removed
%   PACKROWS rewrite
%   RIDGEWALK modified to return "ascending" ridges only
%   _______________________________________________________________________
% 
%   Changes new in version 0.82 
%
%   JLAB_LICENSE slightly altered
%   Minor change to test report convention
%   Missing tests added
%   New functions: PF_PARAMS and PF_EXTRACT
%   New functions: TRACKFILL and TRACK2GRID
%   New function: SLIDETRANS
%   New function: QUADFORM
%   New functions: ZE2DIST and ZE2INC
%   New functions: SPHEREDIST and SPHEREPROJ 
%   New function: RADEARTH
%   New functions: ZEAZ2XY and XY2ZEAZ
%   ELLCONV modified to match new notation, added tests
%   MORSEWAVE changed to specify frequency exactly
%   MORLWAVE changed to exact zero-mean formulation
%   MORSEWAVE and MORLWAVE output argument change 
%   RIDGEWALK output argument change 
%   RIDGEINTERP input and output argument change
%   VINDEXINTO bugfix and testing, led to VSHIFT bugfix 
%   ISSCALAR renamed to ISSCAL to avoid naming conflict
%   COL2MAT bugfix for length of key output matrix
%   VFILT setting filter to unit energy feature removed
%   _______________________________________________________________________ 
%
%   Changes new in version 0.81
%
%   RIDGEINTERP bugfix to have NANs same in frequency
%   LATLON2XY and LATLON2UV have complex NANs for one output argument
%   CATSTRUCT modified to use complex NANs for missing complex data
%   MSPEC functionality split between MSPEC and new function MTRANS
%   STICKVECT bugfix
%   Missing functions included: TIDEFREQ and SPECDIAG
%   _______________________________________________________________________
% 
%   Changes new in version 0.8 
%
%   New functions: XY2LATLON, SPECDIAG, RIDGEINTERP, TIDEFREQ
%   New functions: ECCONV, ELLCONV, ELLDIFF, ELLVEL, ELLRAD, MAKELLIPSE
%   WAVETRANS bugfix for multiple wavelets 
%   VDIFF modified for multiple input arguments
%   MAKE additional input format added 
%   VSHIFT added selection functionality 
%   ELLIPSEPLOT bugfix
%   CIRC2ELL input argument change 
%   WAVESPECPLOT chaged to allow complex-valued trasnform matrices
%   WCONVERT renamed TRANSCONV and syntax changed 
%   COMMENTLINES bugfix for directory arguments
%   RIDGEWALK and RIDGEINTERP modified to detect negative rotary transform
%   LATLON2XY bugfix for larger than 180 degree jumps; test code
%   RUNTESTS_JLAB modified to test modules separately
%   NUMSLABS moved to JARRAY
%   ROT changed to handle special cases of n*pi/2 exactly; test added
%   FILLBAD rewrite, also changed to handle complex-valued data
%   VDIFF changed to differentiate a specified dimension
%   MORSEWAVE changed to have phase definition of frequency
%   VZEROS changed to support NAN output
%   XOFFSET and YOFFSET changed to offset groups of lines
%   UVPLOT streamlined, hold off by default
%   RIDGEWALK argument convention change
%   RIDGEWALK nasty bugfix to form_ridge_chains
%   VFILT bugfix for filter of zeros and ones
%   MAKEFIGS all modified to print to current directory
%   PACKROWS and PACKCOLS changing labels feature removed
%   ELLIPSEPLOT major axis drawing added
%   LATLON2UV bug fixed for vector day input
%   POLPARAM modified for more general input formats
%   MATMULT acceleration and test
%   POLPARAM bugfix for noise causing imaginary determinant 
%   VSUM and VSWAP modified to support possible NAN+i*NAN
%   SLEPWAVE bugfix for 'complex' flag
%   MSPEC rewrite, bugfix for odd length clockwise, test suite
%   _______________________________________________________________________
%
%   Changes new in version 0.72
%
%   New functions: VEMTPY, VCELLCAT, EDGEPOINTS
%   New function:  MODMAXPEAKS
%   New function:  NONNAN
%   BLOCKLEN fixed incorrect output length  
%   ELLIPSEPLOT bug fixed
%   FLUSHLEFT and FLUSHRIGHT re-written
%   LATLON2CART and LATLON2CV renamed LATLON2XY and  LATLON2UV
%        Argument handling improvements to both
%   LATLON2XY swapped order of output arguments
%   MODMAX modified to prevent long traverses 
%   MORSECFUN and MORSEAREA modified for mixed matrix / scalar arguments
%   MORSEWAVE modified to ensure centering of wavelet
%   NORMFORM bug for infinite theta fixed
%   PDFPROPS modified to output skewness and kurtosis
%   RIDGEWALK modified to handle possible absence of ridges
%   RIDGEWALK output argument improvements
%   RIDGEWALK input argument NS removed
%   RIDGEWALK modified to output matrices
%   RIDGERECON modified to handle multivariate datasets
%   RIDGERECON modified to handle ridges of a complex-valued time series
%   SLEPTAP fixed non-unit energy for interpolated tapers
%   TWODHIST fixed counting bug and general improvements
%   WAVETRANS modified for better centering
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002--2011 J.M. Lilly --- type 'help jlab_license' for details

help jlab_changes
