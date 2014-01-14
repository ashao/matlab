%
% JCELL  Tools for operating on cell arrays of numerical arrays.
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, *JCELL*, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Cell / column vector conversions
%   col2cell   - Converts 'column-appended' data into cell arrays of column vectors.
%   cell2col   - Converts cell arrays of column vectors into 'column-appended' data.
%
% Identity, sizes, and extrema
%   cellaresame - Tests whether two cell arrays of arrays are the same.
%   cellength  - Length of each element in a cell array.
%   cellsize   - Size of each element in a cell array along specified dimension.
%   cellmax    - Maximum of each element in a cell array.
%   cellmin    - Minimum of each element in a cell array.
%
% Real, imag, abs, and angle 
%   cellreal   - Real part of each element in a cell array.
%   cellimag   - Imaginary part of each element in a cell array.
%   cellabs    - Absolute value of each element in a cell array.
%   cellangle  - Complex argument (angle) of each element in a cell array.
%   cellconj   - Complex conjugate of each element in a cell array.
%
% Plotting
%   cellplot   - Plot elements of a cell array of numeric arrays.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2011 J.M. Lilly --- type 'help jlab_license' for details      


help Jcell
if 0
           %Cell / column vector conversions
             col2cell   %- Converts 'column%-appended' data into cell arrays of column vectors.
             cell2col   %- Converts cell arrays of column vectors into 'column%-appended' data.
          
           %Identity, sizes, and extrema
             cellaresame %- Tests whether two cell arrays of arrays are the same.
             cellength  %- Length of each element in a cell array.
             cellsize   %- Size of each element in a cell array along specified dimension.
             cellmax    %- Maximum of each element in a cell array.
             cellmin    %- Minimum of each element in a cell array.
          
           %Real, imag, abs, and angle 
             cellreal   %- Real part of each element in a cell array.
             cellimag   %- Imaginary part of each element in a cell array.
             cellabs    %- Absolute value of each element in a cell array.
             cellangle  %- Complex argument (angle) of each element in a cell array.
             cellconj   %- Complex conjugate of each element in a cell array.
          
           %Plotting
             cellplot   %- Plot elements of a cell array of numeric arrays. 
end
