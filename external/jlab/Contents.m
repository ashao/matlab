%
%   JLAB        - Matlab freeware for data analysis 
%   Author      - J. M. Lilly
%   Requires    - Matlab Version 7.10 or higher
%   _________________________________________________________________
%
%   See also: 
%   *CONTENTS*, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
%   See "www.jmlilly.net/jmlsoft.html" for a brief introduction and for
%   the most recent version.
%
%   Send comments, questions, and bug reports to "eponym@jmlilly.net".
%
%   This is copyrighted software and is distributed according to certain
%   very reasonable terms.  See the jlab_license for details.
% 
%   JLAB is organized into a number of modules:
%
%     General purpose
%     -----------------------------------------------------------------
%     Jarray     - Vector, matrix, and N-D array tools.
%     Jmath      - Mathematical aliases and basic functions.
%     Jpoly      - Special polynomials, matrices, and functions.
%     Jgraph     - Fine-turning and customizing figures.
%     Jstrings   - Strings, files, and variables.
%     Jstats     - Statistical tools and probability distributions. 
%     Jsignal    - Signal processing, wavelet and spectral analysis.
%     Jellipse   - Elliptical (bivariate) time series analysis.
%     Jcell      - Tools for operating on cell arrays of numerical arrays.
%     Vtools     - Operations on multiple data arrays simultaneously.
%
%     Special purpose
%     -----------------------------------------------------------------
%     Joceans    - Oceanography-specific functions.
%     Jsphere    - Spherical geometry and derivatives.
%     Jsatfun    - Satellite data treatment and design.
%     Jtriads    - Gravity-capillary triad interaction functions. 
%     Jpapers    - Figures from papers by J. M. Lilly, i.e. myself.
%
%   Type 'help jarray' for contents of the JARRAY module, etc.
%
%   A good starting point is:
%
%    jlab_highlights - The best of JLAB.  You really want to read this.
%
%   Some other important general functions:
%
%     jlab_license   - License statement and permissions for JLAB package.
%     jlab_version   - JLAB version number.
%     jlab_changes   - Changes to JLAB in each release.
%     jlab_runtests  - Run test suite for JLAB package.
%     jlab_settings  - Specifies settings for customizable JLAB properties.
%     jlab_matfiles  - Description of included and optional redistributed data.
%     jlab_acknowledgement  - Sources of support for developing this software.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J. M. Lilly --- type 'help jlab_license' for details        

%   Administrative / hidden
%   jlab_todo     - List of things to do for JLAB, in jlocal.
%   cd_figures - Change directory to figure repository.

if 0
    Jarray;
    Jmath;
    Jpoly;
    Jgraph;
    Jstrings;
    Jstats;
    Jsignal;
    Jellipse;
    Jcell;
    Vtools;
    Joceans;
    Jsphere;
    Jsatfun;
    Jtriads;
    Jpapers;
    jlab_license;
    jlab_version;
    jlab_changes;
    jlab_runtests;
    jlab_settings;
    jlab_matfiles;
    jlab_acknowledgement;
end
