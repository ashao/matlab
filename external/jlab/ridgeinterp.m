function[varargout]=ridgeinterp(varargin)
% RIDGEINTERP  Interpolate quantity values onto ridge locations.
%
%   XI=RIDGEINTERP(W,FS,RQ,IR,JR,X) where W is a wavelet transform at
%   at *radian* frequencies FS, RQ is a ``ridge quantity'', and IR and
%   JR are the time and scale indices of ridges, interpolates quantity 
%   X along the ridge locations and returns the result in XI.
%
%   XI is the same size as IR and JR.  RQ and X are the same size as W.
%
%   W is output by WAVETRAN, RQ is output by RIDGEQUANTITY, and IR and 
%   JR are output by RIDGEWALK.  
%
%   While the transform W is specifed only at discrete frequencies,
%   RIDGEINTERP interpolates transform values between these levels to 
%   find a more precise value of the transform along the ridges than 
%   simply looking up the values of W at rows IR and columns JR.  
%
%   [XI1,XI2,...,XIN]=RIDGEINTERP(W,FS,RQ,IR,JR,X1,X2,...,XN) also 
%   interpolates the quantities X1,X2,...,XN, all the same size as W. 
%
%   RIDGEINTERP uses fast quadratic interpolation via QUADINTERP.  In
%   rare cases, quadratic interpolation fails for an individual point and 
%   therefore linear interpolation is used instead.
%   __________________________________________________________________
%
%   RIDGEINTERP is a low-level function called by RIDGEWALK.
%
%   See also RIDGEWALK, RIDGEQUANTITY, QUADINTERP.
%
%   Usage:  xi=ridgeinterp(w,fs,rq,ir,jr,x);
%           [xi1,xi2,xi3]=ridgeinterp(w,fs,rq,ir,jr,x1,x2,x3);
%           [xi1,xi2,xi3]=ridgeinterp(dt,w,fs,rq,ir,jr,x1,x2,x3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2011 J.M. Lilly --- type 'help jlab_license' for details    


wo=varargin{1};
fs=varargin{2};
rq=varargin{3};
ir=varargin{4};
jr=varargin{5};

args=varargin(6:end);

index=find(~isnan(ir));
if ~isempty(index)
    varargout=ridgeinterp1_quadratic(index,ir,jr,fs,wo,rq,args);
else
    for i=1:length(args)
        varargout{i}=[];
    end
end
              
function[outargs]=ridgeinterp1_quadratic(index,ir,jr,fs,w,rq,args)

sizeir=size(ir);
vcolon(ir,jr);
vindex(ir,jr,index,1);

indexr=nonnan(sub2ind(size(rq),ir,jr));

%Ridge quantity along the ridges, and at one scale up and down

di=size(w,1);
dr=rq(indexr);
drp=rq(indexr+di);
drn=rq(indexr-di);

[xmin,jre]=quadinterp(jr-1,jr,jr+1,abs(drn).^2,abs(dr).^2,abs(drp).^2);

%Rare complete failure of quadratic interpolation associated 
%with the ridge quantity changing sign between jrp and jrn, yet
%not having the ridge quantity at jr being between these two values
        
bool=~( (jr+1>jre) & (jre> jr-1)); 
jre(bool)=lininterp(jr(bool)-1,jr(bool)+1,drn(bool),drp(bool));


for i=1:length(args)
   x=args{i};
   %size(x)
   xr=nan*zeros(sizeir(1),sizeir(2),size(w,3));
   for k=1:size(x,3);
        xk=x(:,:,k);
        
        xro=xk(indexr);
        xrp=xk(indexr+di);
        xrn=xk(indexr-di);
        
        xrk=quadinterp(jr-1,jr,jr+1,xrn,xro,xrp,jre);         
        
        %Use linear interpolation where quadratic fails
        xrk(bool)=lininterp(jr(bool)-1,jr(bool)+1,xrn(bool),xrp(bool),jre(bool)); 
        
       %figure,plot([jr-1 jr+1 jre]) 
           
       %figure,plot(([drn dr drp])) 
       
       %figure,plot(real([xrn xrp xrk]))

        %length(find(bool))
        xr(index,k)=xrk;
   end
   %Linearly interpolate between the approximate ridge and the bracketing curve 
   outargs{i}=xr;
end