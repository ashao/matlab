function[rq,om]=ridgequantity(w,fs,str)
%RIDGEQUANTITY  The ``ridge quantity'' associated with a wavelet transform.
%
%   RQ=RIDGEQUANTITY(W,FS,STR) returns the ridge quantity RQ associated 
%   with the analytic wavelet transform W computed at frequencies FS. RQ is
%   the same size as W.
%
%   STR determines the type of ridge, with STR='amplitdue' for amplitude 
%   ridges and STR='phase' for phase ridges.  
%
%   [RQ,FW]=RIDGEQUANTITY(W,FS,STR) also returns the transform frequency 
%   matrix FW.  FW has the radian instantaneous frequency at each scale and
%   is the same size as W.
% 
%   For details see
%
%       Lilly and Olhede (2010).  On the analytic wavelet transform.  
%           IEEE Trans. Info. Theory.
%   _____________________________________________________________________
%
%   Joint ridges
%
%   IF W has more than one element along its third dimension, then RQ is
%   the joint ridge quantity and F is the joint transform frequency.
%
%   RQ and FW then have the same number of rows and columns as W but only
%   one element along the third dimension.
%
%   For details see
%
%       Lilly and Olhede (2009).  Wavelet ridge estimation of jointly
%       modulated multivariate oscillations. Asilomar Conference on
%       Signals, Systems, and Computers.
%   _____________________________________________________________________
%
%   RIDGEQUANTITY is a low-level function called by RIDGEWALK.
%
%   See also RIDGEWALK, RIDGEINTERP.
%
%   Usage: rq=ridgequantity(w,fs,str);
%          [rq,fw]=ridgequantity(w,fs,str);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details
 
if isstr(w)
    error('No test for RIDGEQUANTITY.')
end
om=instfreq(w);
a=abs(w);    

if size(w,3)~=1
    om=instfreq(w,3);
    a=sqrt(vmean(a.^2,3));
end

fsmat=vrep(fs(:)',size(w,1),1);

if strcmp(str(1:3),'pha')
    rq=om-fsmat;
    %rq=(om-fsmat)./fsmat;
elseif strcmp(str(1:3),'amp')  
    %This is da/ds
    rq=fsmat.^2.*frac(vdiff(a,2),vdiff(fsmat,2)).*frac(1,a);
%elseif strcmp(str(1:3),'gro')
%    dads=frac(fsmat.^2,2*pi).*frac(vdiff(a,2),vdiff(fsmat,2)).*frac(1,a);
%    rq=(om-fsmat)+sqrt(-1)*dads;
% elseif strcmp(str(1:3),'gro')  
%  %    rq=vdiff(sqrt(squared(vdiff(a,2).*frac(1,a))+squared(frac(1,a.*fsmat).*vdiff(a,1))),2);
%    %  rq=sqrt(squared(vdiff(a,2).*frac(1,a))+squared(frac(1,a.*fsmat).*vdiff(a,1)));
%      %rq=frac(1,fsmat).*sqrt(squared(vdiff(a,1).*frac(1,a))+squared(frac(fsmat.^2,a).*frac(vdiff(a,2),vdiff(fsmat,2))));
%      xx=sqrt(squared(vdiff(a,1).*frac(1,a))+squared(frac(fsmat.^2,a).*frac(vdiff(a,2),vdiff(fsmat,2))));
%      rq=-fsmat.^2.*frac(vdiff(xx,2),vdiff(fsmat,2));
%      %rq=xx;
end