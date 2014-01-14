% function [amp,phase,frac,offset,yy]=annual_harms(y,t,nharm,cutoff)
%
% given times t (days) and data in y:
%	computes annual cosine harmonics (n=1:nharm)
%         f(t) = a * cos(2*pi*n*t/L + phi) + b 
%
%  if fraction of variance is > cutoff/100, removes harmonic
%		otherwise, sets amp, etc to zero   
%		in any case, removes offset for first harmonic
%
% Note: cutoff value should be given in PERCENT
%  returns positive amplitude,phase,frac(tion) of variance, and offset
%	for each harmonic
% yy is the residual time series
%
function [yy,amp,phase,frac,offset]=annual_harms(y,t,nharm,cutoff,L)
y=y(:);t=t(:);nt=length(t);
% L=365.25;
%
yy=y;
wt=ones(nt,1);
    for i=1:nharm
    Li=L/i;
    [amp(i),phase(i),offset(i),yi,frac(i)]=fit_cosine(t,yy,wt,Li);
% make amplitude positive
       if amp(i)<0;
       amp(i)=-amp(i);
       phase(i)=phase(i)+pi;
          if phase(i)>pi;phase(i)=phase(i)-2*pi;end
       end
       if frac(i) > cutoff/100
       yy=yi;
       else
       amp(i)=0;
       phase(i)=0;
% remove first offset if first harmonic is not significant
          if i == 1
          yy=y-offset(i)*ones(nt,1);
          end
       end
    end

