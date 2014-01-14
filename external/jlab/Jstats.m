%
% JSTATS   Statistical, distribution, and mapping tools.
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, *JSTATS*, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%   
% Mapping scattered data
%   polysmooth - Smoothing scattered 2D data with local polynomial fitting.
%
% Statistics of data as a function of two independent variables.
%   twodhist   - Two-dimensional histogram.
%   twodmed    - Median value of a function of two variables.
%   twodstats  - Mean, variance, and covariance of functions of two variables.
%
% Sorting etc.
%   bindata    - Rapidly sort data into adjacent bins.
%   twodsort   - Distances from data points to nearby grid points.
%   spheresort - Sorted great circle distances to nearby points on a sphere.
%
% Probability density functions.
%   simplepdf  - Gaussian, uniform, Cauchy, and exponential pdfs.
%
% Moments and cumulants
%   mom2cum    - Convert moments to cumulants.
%   cum2mom    - Convert cumulants to moments.
%
% Operations on probability density functions.
%   conflimit  - Computes confidence limits for a probability density function.
%   pdfprops   - Mean and variance associated with a probability distribution.
%   pdfadd     - Probability distribution from adding two random variables.
%   pdfmult    - Probability distribution from multiplying two random variables
%   pdfdivide  - Probability distribution from dividing two random variables.
%   pdfinv     - Probability distribution of the inverse of a random variable.
%   pdfchain   - The "chain rule" for probabilty density functions.
%   pdfconv    - Convolution of a probability distribution with itself.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details      

help Jstats

if 0              
           %Mapping scattered data
             polysmooth %- Smoothing scattered 2D data with local polynomial fitting.
          
           %Statistics of data as a function of two independent variables.
             twodhist   %- Two%-dimensional histogram.
             twodmed    %- Median value of a function of two variables.
             twodstats  %- Mean, variance, and covariance of functions of two variables.
          
           %Sorting etc.
             bindata    %- Rapidly sort data into adjacent bins.
             twodsort   %- Distances from data points to nearby grid points.
             spheresort %- Sorted great circle distances to nearby points on a sphere.
          
           %Probability density functions.
             simplepdf  %- Gaussian, uniform, Cauchy, and exponential pdfs.
          
           %Moments and cumulants
             mom2cum    %- Convert moments to cumulants.
             cum2mom    %- Convert cumulants to moments.
          
           %Operations on probability density functions.
             conflimit  %- Computes confidence limits for a probability density function.
             pdfprops   %- Mean and variance associated with a probability distribution.
             pdfadd     %- Probability distribution from adding two random variables.
             pdfmult    %- Probability distribution from multiplying two random variables
             pdfdivide  %- Probability distribution from dividing two random variables.
             pdfinv     %- Probability distribution of the inverse of a random variable.
             pdfchain   %- The "chain rule" for probabilty density functions.
             pdfconv    %- Convolution of a probability distribution with itself.
end