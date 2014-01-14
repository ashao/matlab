%
% JARRAY   Vector, matrix, and N-D array tools.
%
%   See also:
%   CONTENTS, *JARRAY*, JMATH, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Vector and matrix multiplication
%   oprod         - Outer product of two column vectors. 
%   osum          - Outer sum of two column vectors.
%   vectmult      - Matrix multiplication for arrays of N-vectors.  
%   offset        - Offsets matrix rows or columns.
%   matmult       - Matrix multiplication for arrays of matrices. 
%   matinv        - Fast inversion of arrays of small matrices.
%
% Array transformation and reduction
%   row2col       - Given any vector, returns it as a column vector.
%   col2row       - Given any vector, returns it as a row vector.
%   nonnan        - Return all non-NAN elements of an array.
%  
% Arrays of blocks of non-uniform length
%   blocknum      - Numbers the contiguous blocks of an array.
%   blocklen      - Counts the lengths of 'blocks' in an array.
%
% Multidimensional (N-D) array tools
%   nd            - Number of the last nonsingleton dimension.
%   nnsd          - Number of nonsingleton dimensions.
%   ndtrans       - Generalized transpose of a potentially multidimensional vector.
%   ndindex       - Indexes a multidimensional array along a specified dimension
%
% Array tests
%   aresame       - Tests whether two N-D arrays are the same.
%   iseven        - Tests whether the elements of an array are even.
%   isodd         - Tests whether the elements of an array are odd.
%   ismat         - Tests whether the argument is a 2-D matrix; false for scalars.
%   issing        - Tests whether the argument is a singleton array.
%   issquare      - Tests whether the argument is a square matrix.
%   jisrow        - Tests whether the argument is a row vector.
%   jiscol        - Tests whether the argument is a column vector.
%   islargest     - True for the largest-magnitude element at each index location.
%   turningpoint  - True for turning points, i.e. local extrema, along rows.
%
% Array operation aliases -- min, max, etc.
%   allall        - ALLALL(X)=ALL(X(:))
%   maxmax        - MAXMAX(X)=MAX(X(:))
%   minmin        - MINMIN(X)=MIN(X(:))
%   anyany        - ANYANY(X)=ANY(X(:))
%   sumsum        - SUMSUM(X)=SUM(X(:))
%   range         - RANGE(X)=[MIN(X(:)) MAX(X(:))]
%
% Indexing (subscripting) and sorting
%   crossings     - Find values of an array crossing a specified threshold.
%   twodsort      - Distances from data points to nearby grid points.
%
% Interpolation
%   jinterp       - Matrix-matrix 1-D interpolation.
%   doublen       - Interpolates a time series to double its length.
%   sinterp       - Spline-interpolates a column vector to a new length.
%   fillbad       - Linearly interpolate over bad data points. 
%   lininterp     - Fast linear interpolation for arbitrary-sized arrays.
%   quadinterp    - Fast quadratic interpolation for arbitrary-sized arrays.
%   cubeinterp    - Fast cubic interpolation for arbitrary-sized arrays.
%
% Smoothing
%   smartsmooth   - Fast light smoothing for large matrices.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details        


help Jarray

if 0
              
           %Vector and matrix multiplication
             oprod         %-     Outer product of two column vectors. 
             osum          %-     Outer sum of two column vectors.
             vectmult      %-     Matrix multiplication for arrays of N%-    vectors.  
             offset        %-     Offsets matrix rows or columns.
             matmult       %-     Matrix multiplication for arrays of matrices. 
             matinv        %-     Fast inversion of arrays of small matrices.
          
          % Array transformation and reduction
             row2col       %-     Given any vector, returns it as a column vector.
             col2row       %-     Given any vector, returns it as a row vector.
             nonnan        %-     Return all non%-    NAN elements of an array.
             finite        %-     Return all finite elements of an array.
            
           %Arrays of blocks of non%-    uniform length
             blocknum      %-     Numbers the contiguous blocks of an array.
             blocklen      %-     Counts the lengths of 'blocks' in an array.
          
          % Multidimensional (N%-    D) array tools
             nd            %-     Number of the last nonsingleton dimension.
             nnsd          %-     Number of nonsingleton dimensions.
             ndtrans       %-     Generalized transpose of a potentially multidimensional vector.
             ndindex       %-     Indexes a multidimensional array along a specified dimension
          
           %Array tests
             aresame       %-     Tests whether two N%-    D arrays are the same.
             iseven        %-     Tests whether the elements of an array are even.
             isodd         %-     Tests whether the elements of an array are odd.
             ismat         %-     Tests whether the argument is a 2%-    D matrix; false for scalars.
             issing        %-     Tests whether the argument is a singleton array.
             issquare      %-     Tests whether the argument is a square matrix.
             jisrow        %-     Tests whether the argument is a row vector.
             jiscol        %-     Tests whether the argument is a column vector.
             islargest     %-     True for the largest%-    magnitude element at each index location.
             turningpoint  %-     True for turning points, i.e. local extrema, along rows.
          
           %Array operation aliases %-    %-     min, max, etc.
             allall        %-     ALLALL(X)=ALL(X(:))
             maxmax        %-     MAXMAX(X)=MAX(X(:))
             minmin        %-     MINMIN(X)=MIN(X(:))
             anyany        %-     ANYANY(X)=ANY(X(:))
             sumsum        %-     SUMSUM(X)=SUM(X(:))
             range         %-     RANGE(X)=[MIN(X(:)) MAX(X(:))]
          
           %Indexing (subscripting) and sorting
             crossings     %-     Find values of an array crossing a specified threshold.
             twodsort      %-     Distances from data points to nearby grid points.
          
           %Interpolation
             jinterp       %-     Matrix%-    matrix 1%-    D interpolation.
             doublen       %-     Interpolates a time series to double its length.
             sinterp       %-     Spline%-    interpolates a column vector to a new length.
             fillbad       %-     Linearly interpolate over bad data points. 
             lininterp     %-     Fast linear interpolation for arbitrary%-    sized arrays.
             quadinterp    %-     Fast quadratic interpolation for arbitrary%-    sized arrays.
             cubeinterp    %-     Fast cubic interpolation for arbitrary%-    sized arrays.
          
           %Smoothing
             smartsmooth   %-     Fast light smoothing for large matrices.
end
