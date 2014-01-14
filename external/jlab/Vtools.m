%
% VTOOLS  Operations on multiple data arrays simultaneously.
%
%   See also:
%   CONTENTS, JARRAY, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, *VTOOLS*, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Adding, multiplying, filtering
%   vadd       - Vector-matrix addition without "dimensional" hassle. 
%   vmult      - Vector-matrix multiplication without "dimensional" hassle.
%   vfilt      - Filtering along rows without change in length.
%   vdiff      - Length-preserving first central difference. 
%   vpower     - Raises array to the specified power.
%   
% Dimensions, sums, moments
%   vnd        - Number of dimensions. 
%   vsize      - Returns the sizes of multiple arguments.
%   vsum       - Sum over finite elements along a specified dimension.
%   vmean      - Mean over finite elements along a specified dimension.
%   vstd       - Standard deviation over finite elements along a specfied dimension.
%   vmoment    - Central moment over finite elements along a specfied dimension.    
%   vmedian    - Median over finite elements along a specified dimension.
%
% Indexing, shifting, swapping
%   vindex     - Indexes an N-D array along a specified dimension. 
%   vrep       - Replicates an array along a specified dimension.                   
%   vshift     - Cycles the elements of an array along a specified dimension.       
%   vsqueeze   - Squeezes multiple input arguments simultaneously. 
%   vswap      - Swap one value for another in input arrays.
%   vcolon     - Condenses its arguments, like X(:).                           
%   vcellcat   - Concatenates cell arrays of column vectors.
%
% Initializing
%   vempty     - Initializes multiple variables to empty sets.
%   vzeros     - Initializes multiple variables to arrays of zeros.
%
% Datasets of non-uniform length: column / matrix / cell conversions. 
%   colbreaks     - Insert NANs into discontinuties in a vector.
%   mat2col       - Compress NAN-padded matrix data into long columns.
%   col2mat       - Expands 'column-appended' data into a matrix.
%   cell2col      - Converts cell arrays of column vectors into 'column-appended' data.
%   col2cell      - Converts 'column-appended' data into cell arrays of column vectors.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details        

help Vtools

if 0          
           %Adding, multiplying, filtering
             vadd       %-        Vector-matrix addition without "dimensional" hassle. 
             vmult      %-        Vector-matrix multiplication without "dimensional" hassle.
             vfilt      %-        Filtering along rows without change in length.
             vdiff      %-        Length-preserving first central difference. 
             vpower     %-        Raises array to the specified power.
             
           %Dimensions, sums, moments
             vnd        %-        Number of dimensions. 
             vsize      %-        Returns the sizes of multiple arguments.
             vsum       %-        Sum over finite elements along a specified dimension.
             vmean      %-        Mean over finite elements along a specified dimension.
             vstd       %-        Standard deviation over finite elements along a specfied dimension.
             vmoment    %-        Central moment over finite elements along a specfied dimension.    
             vmedian    %-        Median over finite elements along a specified dimension.
          
           %Indexing, shifting, swapping
             vindex     %-        Indexes an N-D array along a specified dimension. 
             vrep       %-        Replicates an array along a specified dimension.                   
             vshift     %-        Cycles the elements of an array along a specified dimension.       
             vsqueeze   %-        Squeezes multiple input arguments simultaneously. 
             vswap      %-        Swap one value for another in input arrays.
             vcolon     %-        Condenses its arguments, like X(:).                           
             vcellcat   %-        Concatenates cell arrays of column vectors.
          
          % Initializing
             vempty     %-        Initializes multiple variables to empty sets.
             vzeros     %-        Initializes multiple variables to arrays of zeros.
          
           %Datasets of non-uniform length: column / matrix / cell conversions. 
             colbreaks     %-        Insert NANs into discontinuties in a vector.
             mat2col       %-        Compress NAN-padded matrix data into long columns.
             col2mat       %-        Expands 'column-appended' data into a matrix.
             cell2col      %-        Converts cell arrays of column vectors into 'column-appended' data.
             col2cell      %-        Converts 'column-appended' data into cell arrays of column vectors.
end
             
             
             