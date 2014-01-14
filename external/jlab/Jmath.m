%
% JMATH  Mathematical aliases and basic functions.
%
%   See also:
%   CONTENTS, JARRAY, *JMATH*, JPOLY, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Mathematical aliases -- algebraic shorthand   
%   oneover       - Returns reciprocal of argument.
%   frac          - Makes a fraction: FRAC(A,B)=A./B
%   twopi         - Raise (2 pi) to the specified power.
%   res           - Residual after flooring:  RES(X)=X-FLOOR(X)
%   squared       - Squares its argument:  SQUARED(X)=X.^2
%   cdot          - CDOT(A,B)=REAL(A.*CONJ(B))
%   rot           - Complex-valued rotation:  ROT(X)=EXP(SQRT(-1)*X)
%   choose        - Binomial coefficient: CHOOSE(N,K) = "N -choose- K"
%
% Boolean expressions, sets, and slabs
%   subset        - Extract subset of A given B.  
%   lookup        - Locate elements of one array within another.
%   ismemb        - Tests whether the elements of an array are members of a set.
%   iscompat      - Tests whether an array's size is "compatible" with another's.
%   numslabs      - Counts the number of 'slabs' of one array relative to another.
%
% Other
%   quadform      - Implements the quadratic formula.
%
% Angles, degrees, and radians
%   relog         - RELOG(X)=REAL(LOG(X))
%   imlog         - IMLOG(X)=UNWRAP(IMAG(LOG(X)))
%   jrad2deg      - Converts radians to degrees.
%   jdeg2rad      - Converts degrees to radians.
%   deg360        - Converts degrees to the range [0, 360].
%   deg180        - Converts degrees to the range [-180, 180].
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details      

help Jmath

if 0           
           %Mathematical aliases %-%- algebraic shorthand   
             oneover       %- Returns reciprocal of argument.
             frac          %- Make fraction: FRAC(A,B)=A./B
             twopi         %- Raise (2 pi) to the specified power.
             res           %- Residual after flooring:  RES(X)=X%-FLOOR(X)
             squared       %- Squares its argument:  SQUARED(X)=X.^2
             cdot          %- CDOT(A,B)=REAL(A.*CONJ(B))
             rot           %- Complex%-valued rotation:  ROT(X)=EXP(SQRT(%-1)*X)
             choose        %- Binomial coefficient: CHOOSE(N,K) = "N choose K"
          
           %Boolean expressions, sets, and slabs
             subset        %- Extract subset of A given B.  
             lookup        %- Locate elements of one array within another.
             ismemb        %- Tests whether the elements of an array are members of a set.
             iscompat      %- Tests whether an array's size is "compatible" with another's.
             numslabs      %- Counts the number of 'slabs' of one array relative to another.
          
           %Other
             quadform      %- Implements the quadratic formula.
          
           %Angles, degrees, and radians
             relog         %- RELOG(X)=REAL(LOG(X))
             imlog         %- IMLOG(X)=UNWRAP(IMAG(LOG(X)))
             jrad2deg      %- Converts radians to degrees.
             jdeg2rad      %- Converts degrees to radians.
             deg360        %- Converts degrees to the range [0, 360].
             deg180        %- Converts degrees to the range [%-180, 180]. 
end
