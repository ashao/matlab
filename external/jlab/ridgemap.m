function[varargout]=ridgemap(varargin)
%RIDGEMAP  Maps ridge quantities back onto the time series.
%
%   X=RIDGEMAP(N,XR,IR) where IR is a ridge index and XR is a quantity 
%   along the ridge, maps the values of XR to their correct row locations 
%   in a time series of length N, and returns the result in the array X.
%
%   If IR and XR contain L ridges separated by NANs, as output by
%   RIDGEWALK, then X is N x L with the values of XR from each ridge in a 
%   separate column. Values not specified by the IR are left as NANs. 
%
%   [X1,X2,...,XM]=RIDGEMAP(N,X1R,X2R,...,XPR,IR) also works for any P
%   ridge quantities X1R--XPR.
%
%   The XR and IR quantities are assumed to have been created by RIDGEWALK.
%   __________________________________________________________________
% 
%   Joint ridges
%
%   X=RIDGEMAP(SIZ,XR,IR) works for joint ridges.  In this case the
%   original signal has size SIZ=[M N].  Then XR will have M columns, one 
%   for each time series component, and the output X is N x M x L. 
%
%   [X1,X2,...,XM]=RIDGEMAP(SIZ,X1R,X2R,...,XPR,IR) also works here.
%
%   See RIDGEWALK for more details on joint ridges.
%   __________________________________________________________________
%
%   Collapsing 
%
%   X=RIDGEMAP(...'collapse') combines values from all the ridges using
%   a power-weighted mean.  Then X is N x 1, or N x M for joint ridges.
%   __________________________________________________________________
% 
%   Ridge multiplicity
%
%   [...,MULT]=RIDGEMAP returns the ridge multiplicity MULT after all the
%   expected output quantities.  MULT is a column vector with size N x 1.
% 
%   The ridge multiplicity is the number of ridges present at each time.     
%   __________________________________________________________________
%
%   See also RIDGEWALK, RIDGEINTERP.
%
%   Usage:   x=ridgemap(N,xr,ir);
%            [x,f]=ridgemap(N,xr,fr,ir);
%            [x,f]=ridgemap(N,xr,fr,ir,'collapse');
%            [x,mult]=ridgemap(N,xr,ir);
%            [x,f,mult]=ridgemap(N,xr,fr,ir);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009--2011 J.M. Lilly --- type 'help jlab_license' for details

if strcmp(varargin{1}, '--t')
    ridgemap_test,return
end

if isstr(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='all';
end

N=varargin{1}(1);
if length(varargin{1})==1
    M=1;
else
    M=varargin{1}(2);
end

varargin=varargin(2:end);
ir=varargin{end};
varargin=varargin(1:end-1);

L=length(find(isnan(ir)));
if L==0
    L=1;
end

%Initialize output variables
for i=1:length(varargin)
    if ~isreal(varargin{i})
       varargout{i}=(nan+sqrt(-1)*nan)*ones(N,M,L);
    else
       varargout{i}=nan*ones(N,M,L);  
    end
end

if ~isempty(ir)
    for i=1:length(varargin)
        temp=[];
        for k=1:M
            [irtemp,temp(:,k,:)]=col2mat(ir,varargin{i}(:,k));
        end
        varargin{i}=temp;
    end
    ir=col2mat(ir);

    for i=1:length(varargin)
        for j=1:M
            for k=1:size(ir,2)
                 varargout{i}(nonnan(ir(:,k)),j,k)=varargin{i}(find(nonnan(ir(:,k))),j,k);
            end
        end
    end
end


%Calculate multiplicity
%mult=~isnan(vswap(varargout{1}(:,1,:),0,nan));
mult=vsum(~isnan(varargout{1}(:,1,:)),3);
varargout{end+1}=mult;

if strfind(str,'col')        
    for i=length(varargin):-1:1
        if i==1
            varargout{i}=vsum(varargout{i},3);
        else
            varargout{i}=powermean(varargout{i},varargout{1},3);
        end
    end
end

for i=1:length(varargin)        
    varargout{i}=squeeze(varargout{i});
end


function[]=ridgemap_test

load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp,wn]=vectmult(tmat,wx,wy);

%Form ridges of component time series
[ir,jr,wr,fr]=ridgewalk(dt,wn,fs,{0,0,'phase'}); 

[wa,fa,mult]=ridgemap(length(wn),wr,fr,ir);
reporttest('RIDGEMAP has one column per ridge, non-joint ridges',size(fa,2)==size(wa,2)&&size(fa,2)==length(find(isnan(ir))))

[ir,jr,wr,fr]=ridgewalk(dt,wp,wn,fs,{0,0,'amp'});   
[wa,fa,mult]=ridgemap([length(wn) 2],wr,fr,ir);
reporttest('RIDGEMAP has one page per ridge, joint ridges',size(fa,3)==size(wa,3)&&size(fa,3)==length(find(isnan(ir))))
reporttest('RIDGEMAP has one column per component, joint ridges',size(fa,2)==size(wa,2)&&size(fa,2)==size(wr,2))


