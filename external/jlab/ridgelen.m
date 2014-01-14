function[lr]=ridgelen(varargin)
%RIDGELEN  Wavelet ridge length expressed as number of full cycles.
%
%   LEN=RIDGELEN(IR,JR,FR) determines the length of the ridges given by
%   IR, JR, and FR, with LEN expressed in number of cycles completed along
%   the ridge.
%
%   IR and JR are arrays of ridge indices, and FR is the frequency of the 
%   wavelet transform value along the ridge, as output by RIDGEWALK. 
%
%   LEN is a row vector with the same size as the number of columns as the 
%   input arrays.
%
%   Usage: lr=ridgelen(ir,jr,fr);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2011 J.M. Lilly --- type 'help jlab_license' for details

%   RIDGELEN is a low-level function called by RIDGEWALK.

args=varargin;

%In call from within RIDGEWALK, ID is already known and the call is 
%lr=ridgelen(id,ir,jr,fr);
if nargin==4
    id=args{1};
    args=args(2:end);
elseif nargin==3
    id=[];
end

ir=args{1};
jr=args{2};
fr=args{3};
    
if iscell(ir)
    for k=1:length(ir)
        if isempty(id)
            idk=[];
        else
            idk=id{k};
        end
        lr{k}=ridgelenloop(idk,ir{k},jr{k},fr{k});   
    end
else
    lr=ridgelenloop(id,ir,jr,fr);   
end

function[lr]=ridgelenloop(id,ir,jr,fr)    
if isempty(id)
    id=cumsum(isnan(ir),1);
end
index=~isnan(ir);
lr=nan*ir;
if ~isempty(index)
    lr(index)=ridgelen1(id(index),ir(index),jr(index),fr(index));
end
    
function[lr]=ridgelen1(id,ir,jr,fr) 
fr=fr./(2*pi);  %Convert to cyclic frequency
 
[num,a,b]=blocknum(id);
%deal with start and end nans
%npoints=~isnan(fr);
vswap(fr,nan,0);
%angr=cumsum(fr.*dt,1);
ar=cumsum(fr,1);  %Ridge age
lena=abs(ar(b)-ar(a));

len1=zeros(size(id));
len1(a)=lena;
len1=cumsum(len1);

len2=zeros(size(id));
len2(a)=[0;lena(1:end-1)];
len2=cumsum(len2);

lr=len1-len2;  %maximum age of each ridge