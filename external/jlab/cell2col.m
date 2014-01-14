function[varargout]=cell2col(varargin)
%CELL2COL  Converts cell arrays of column vectors into 'column-appended' data.
%
%   COL=CELL2COL(X) converts X, a cell array of N column vectors, into a 
%   single column vector COL with N blocks of data, each ending with a NAN.
%
%   X should contain no NANs.  Missing or bad data may be marked by INFs.
%   __________________________________________________________________
%
%   Multiple input /output arguments
%
%   [C1,C2,...,CN]=CELL2COL(X1,X2,...,XN) also works.  
%   __________________________________________________________________
% 
%   Variable overwriting
%
%   CELL2COL(X1,X2,...,XN) with no output arguments overwrites the input
%   variables.
%   __________________________________________________________________
%
%   Invertibility
%  
%   CELL2COL is inverted by COL2CELL, provided 
%
%       (i)   The cell arrays XN are all the same size and
%       (ii)  The first input argument X1 contains no NANs.
%   __________________________________________________________________
%
%   See also COL2CELL, COL2MAT, MAT2COL, COLBREAKS, VCELLCAT.
%
%   Usage: col=cell2col(x);
%          [c1,c2,...,cN]=cell2col(x1,x2,...,xN);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2009 J.M. Lilly --- type 'help jlab_license' for details

if strcmp(varargin{1},'--t')  
   cell2col_test;return
end
bnan=0;
for i=1:nargin
    ray=varargin{i};
    for j=1:length(ray)
        if i==1
             if anyany(isnan(ray{j}))
                bnan=1;
             end
        end
        ray{j}=[ray{j};nan];
    end
    varargout{i}=vcellcat(ray);
end
if bnan
    disp('Warning: CELL2COL found interior NANs; better to use INFs for missing data.')
end        
eval(to_overwrite(nargin));

function[]=cell2col_test
%Tests for cell2col are contained in col2cell