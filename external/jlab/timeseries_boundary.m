function[y]=timeseries_boundary(x,str,bdetrend)
%TIMESERIES_BOUNDARY  Apply boundary conditions to data before transform.
%
%   TIMESERIES_BOUNDARY applies periodic, zero-padded, or mirror boundary
%   conditions to a time series.  See ANATRANS or WAVETRANS.
%
%   This is a low-level function, not meant to be user called by users.
%
%   Usage: y=timeseries_boundary(x,str,detrendstr);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 
%Prepare data by applying boundary condition

for i=1:size(x,2)
  ai=find(~isnan(x(:,i)),1,'first');
  bi=find(~isnan(x(:,i)),1,'last');
  if isempty(ai) && isempty(bi)
    disp(['Warning: Data column ', int2str(i), ' contains no finite values.'])
    a(i)=1;
    b(i)=size(x,1);
  elseif any(~isfinite(x(ai:bi,i)))
    error(['Data contains interior NANs or INFs in column ', int2str(i), '.'])
  else
    a(i)=ai;
    b(i)=bi;
  end
end

M=size(x,1);
N=size(x,2);

if strcmp(str,'periodic')
   y=zeros(M,N);
else
   y=zeros(3*M,N);
end


for i=1:size(x,2)
   index{i}=a(i):b(i);
   indexy{i}=(M+a(i)-length(index{i}):M+2*length(index{i})+a(i)-1);
   xi=x(index{i},i);

   if bdetrend
       xi=detrend(xi);
   end

   if strcmp(str,'zeros')
       y(indexy{i},i)=[0*xi;xi;0*xi];
   elseif strcmp(str,'mirror')
       y(indexy{i},i)=[flipud(xi);xi;flipud(xi)];
   elseif strcmp(str,'periodic')
       y(index{i},i)=xi;
   else 
       error(['Transform option STR = ''',str,''' is not supported.']);
   end   
end
y=vswap(y,nan,0);

