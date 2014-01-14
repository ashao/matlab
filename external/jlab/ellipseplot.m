function[h]=ellipseplot(varargin)
%ELLIPSEPLOT  Plot ellipses.
%
%   ELLIPSEPLOT(K,L,TH,X) plots an ellipse of amplitude K, linearity L, and
%   orientation TH at complex-valued location X.
%
%   Multiple ellipses are plotted if K, L, TH, and X, are matrices of the
%   the same size. 
%
%   ELLIPSEPLOT draws each column of the input arrays with a different 
%   color, cycling through the default line colors.
%
%   H=ELLIPSEPLOT(...) outputs an array of line handles.
%   ______________________________________________________________________
%
%   Aspect ratio
%
%   ELLIPSEPLOT(...,AR) with AR=[XAR YAR] multiplies the X-signal and
%   Y-signal by XAR and YAR, respectively, for plotting purposes. The
%   aspect ratio is set such that circles appear circular.
%
%   It's okay for AR to have three entries, AR=[XAR YAR ZAR], as is
%   output by "get(gca,'dataaspectratio')".  The last entry is ignored.
%
%   AR is optional and defaults to [1 1]. 
%   ______________________________________________________________________
%
%   Skipping ellipses
%
%   Frequently one does not wish to plot all the ellipses.  There are two
%   ways to accomplish this, as documented in this and the next section.
%
%   K, L, TH, and X are column vectors, plot every SKIP-th ellipse with
%
%         ELLIPSEPLOT(K,L,TH,X, ... ,'skip',SKIP).
%
%   This plots ellipses at indicies [SKIP:SKIP:LENGTH(K)-SKIP].
%
%   More generally, if the input field are matrices of the same size, each
%   having N dimensions, then SKIP can be an array with N elements.  Then
%   SKIP(1) indicates the SKIP number for the first matrix dimension, etc.
%   ______________________________________________________________________
%
%   Indexing ellipses
%   
%   To have more precise control over which ellipses are plotted, use
%   
%          ELLIPSEPLOT(K,L,TH,X, ... ,'index',INDEX) 
%
%   which only plots the ellipses at the indicies indicated by INDEX. If
%   INDEX is empty, nothing happens.  
%
%   When K, L, TH, and X are column vectors, INDEX is an array.  
% 
%   More generally, INDEX can be a cell array of N arrays, one for each 
%   dimenion of the input matrices.
%
%   See PERIODINDEX for generating an index that skips every N periods.
%   ______________________________________________________________________
% 
%   Additional options
%
%   The following trailing options can occur in any order, as long as 
%   they are after the numerical arguments listed above.
%
%   ELLIPSEPLOT(K,L,TH,X, ... ,'phase',PHI) optionally draws a small line, 
%   like a clock hand, to indicate the ellipse phase PHI.
%
%   ELLIPSEPLOT(K,L,TH,X, ... ,'axis') alternatively draws the major axis.
%
%   ELLIPSEPLOT(K,L,TH,X, ... ,'npoints',N) plots ellipses with N points
%   along the circumference.  The default value is 32.  Use N=16 or 64 
%   for faster plotting or for smoother, more well-defined ellipses.
%
%   ELLIPSEPLOT(K,L,TH,X, ... ,STY) also works, where STY is a style string
%   in LINESTYLE format, or a cell array of such strings.  See LINESTYLE.    
%   _______________________________________________________________ 
%
%   See also PERIODINDEX.
%
%   Usage: ellipseplot(k,l,th,x)
%          ellipseplot(k,l,th,x,ar)
%          ellipseplot(k,l,th,x,'axis')
%          ellipseplot(k,l,th,x,ar,'phase',phi)
%          ellipseplot(k,l,th,x,ar,'npoints',64)
%          ellipseplot(k,l,th,x,ar,'index',index)
%          ellipseplot(k,l,th,x,ar,'skip',5)
%          ellipseplot(k,l,th,x,'2r--')
%
%   'ellipseplot --f' generates a sample figure.
%   ______________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details        

if strcmp(varargin{1},'--f')
  ellipseplot_fig;return
end

%   ELLIPSEPLOT(K,L,TH,DX) where DX is a scalar uses DX as a complex-valued
%   offset in between ellipses, beginning at 0. 

%/********************************************************
%Sort out input arguments
na=length(varargin);

k=varargin{1};
l=varargin{2};

[a,b]=kl2ab(k,l);

th=varargin{3};
x=zeros(size(a));
phi=nan*zeros(size(a));
baxis=0;
npoints=32;
%sty{1}='b';sty{2}='g';sty{3}='r';sty{4}='c';sty{5}='m';sty{6}='y';sty{7}='k';
sty{1}='b';
str='matlab';


naold=na+1;
index{1}=-1;
skip=[];
while naold~=na
    naold=na;   
    if ischar(varargin{end})
        if  ~strcmp(varargin{end},'axis')&~strcmp(varargin{end},'mat')&~strcmp(varargin{end},'m_p')
            clear sty
            sty{1}=varargin{end};
        elseif strcmp(varargin{end},'axis')
            baxis=1;
        else 
            str=varargin{end};
        end       
        varargin=varargin(1:end-1);
        na=na-1;     
    end
       
    if iscell(varargin{end})
       sty=varargin{end};
       varargin=varargin(1:end-1);
       na=na-1;  
    end
    
    if ischar(varargin{end-1})
      if strcmp(varargin{end-1},'phase')
          phi=varargin{end};
      elseif strcmp(varargin{end-1},'npoints')
          npoints=varargin{end};
      elseif strcmp(varargin{end-1},'index')
          temp=varargin{end};
          if iscell(temp)
              index=temp;
          else
              index{1}=temp;
          end
      elseif strcmp(varargin{end-1},'skip')
          skip=varargin{end};
      end
      na=na-2;
      varargin=varargin(1:end-2);
    end
end

if na>3
  x=varargin{4};
end


if na>4
  ar=varargin{5};
  if size(ar,2)~=2&&size(ar,2)~=3
    error('The aspect ratio AR must have length equal to 2 or 3.')
  end
  if length(ar)>3,warning('I was expecting AR to be length 2 or 3.'),end
  ar1=ar(1);
  ar2=ar(2);
else
  ar1=1;
  ar2=1;
end
%\********************************************************


if ~isempty(skip)
    for i=1:length(skip)
        index{i}=skip(i):skip(i):size(a,i)-skip(i);
    end
end

h=[];
for i=1:length(index)
    if isempty(index{i})
        disp('ELLIPSEPLOT no points specified, plotting nothing.')
        return  %Yes I am using a return here
    end
end

% N=size(a,1);
% if length(x)==1  &&  N>1
%     x=exp(sqrt(-1)*angle(x)).*conj((0:abs(x):(N-1)*abs(x))');
%     x=oprod(x,ones(size(a(1,:)))');
% end

if index{1}~=-1
    for i=1:length(index)
        vindex(a,b,th,phi,x,index{i},i);
    end
end
clear index
id=osum(zeros(size(a(:,1))),(1:size(a,2))');
vcolon(a,b,th,phi,x,id);
vswap(a,0,nan);
index=find(~isnan(a));

bhold=ishold;
hold on

%Make sure signal will be complex-valued
b(b==0)=1e-10;

storestate=get(gcf,'BackingStore');
set(gcf,'BackingStore','off')
if ~isempty(index)     
    vindex(a,b,th,phi,x,id,index,1);
    h=ellipseplot1(a,b,th,phi,x,ar1,ar2,baxis,npoints,str);
    for i=1:length(sty)
         %I can deal with this by forming a long string like '2b','2b', etc
         hindex=res((id-i)./length(sty))==0;
         if ~isempty(hindex)
            linestyle(h(hindex),sty{i});
         end
    end
end

set(gca,'dataaspectratio',[ar1 ar2 1])
set(gcf,'BackingStore',storestate)

if ~bhold
    hold off
end
set(gca,'box','on')

if nargout==0
  clear h
end

function[h]=ellipseplot1(a,b,th,phi,x,ar1,ar2,baxis,npoints,str)
%z1=(0:.1:2*pi+.1)';

L=size(a,1);
z1=rot(linspace(0,2*pi+0.01,npoints)');
 

if ~allall(isnan(phi))
  z1=[0+sqrt(-1)*0;z1];
end

if baxis
  z1=[-1+sqrt(-1)*0;0+sqrt(-1)*0;z1];
end

x=osum(0*z1,x);
z=osum(z1,0*th);
z=circ2ell(z,a,b,th,phi,ar1,ar2)+x;
if strcmp(str(1:3),'mat')
      h=plot(z);%set(gca,'dataaspectratio',[ar(:)' 1]);
elseif strcmp(str(1:3),'m_p')
      h=m_plot(real(z),imag(z));
end

function[z]=circ2ell(z,a,b,th,phi,ar1,ar2)
%CIRC2ELL  Converts a complex-valued circle into an ellipse.
%
%   ZP=CIRC2ELL(Z,K,L,TH,PHI) where Z is a complex-valued time series,
%   performs a combined phase-lag, stretching, and rotation of Z.  
%  
%   If Z is a circle expressed as a complex-valued time series,
%   e.g. as output by PHASECIRCLE, TH and PHI are angles in radians,
%   and A and B are real-valued weighting factors, CIRC2ELL transforms
%   Z into an ellipse specified by ZP.
%  
%   This ellipse has major axis A and minor axis B, with the major
%   axis oriented at angle TH measured counterclockwise with respect
%   to the positive real axis.  If the phase at temporal midpoint of Z
%   is zero, than the ellipse has phase PHI at the midpoint.
%
%   The input arguments TH, A, B, and PHI may each be scalars or arrays
%   of the same size as Z.
%
%   ZP=CIRC2ELL(... AR) optionally rescales the ellipse by aspect ratio
%   AR for plotting purposes.
%
%   See Lilly (2005) for algorithm details.
%
%   See also ELLIPSEPLOT, PHASECIRCLE
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        
  

rz=real(z);
iz=imag(z);

if isscalar(a)
  a=a+0*z(1,:);
end

if isscalar(b)
  b=b+0*z(1,:);
end

if isscalar(th)
  th=th+0*z(1,:);
end

if isscalar(phi)
  phi=phi+0*z(1,:);
end

phi=vswap(phi,nan,0);

for i=1:size(z,2)
  [rz(:,i),iz(:,i)]=vectmult(jmat(phi(i)),rz(:,i),iz(:,i));
  rz(:,i)=rz(:,i).*a(i);
  iz(:,i)=iz(:,i).*b(i);
  [rz(:,i),iz(:,i)]=vectmult(jmat(th(i)),rz(:,i),iz(:,i));
  if length(ar1)>1   
    rz(:,i)=rz(:,i).*ar1(i);
    iz(:,i)=iz(:,i).*ar2(i);
  else
    rz(:,i)=rz(:,i).*ar1;
    iz(:,i)=iz(:,i).*ar2;
  end
end

z=rz+sqrt(-1).*iz;
%z=vrep(ar1',size(z,1),1).*real(z)+vrep(ar2',size(z,1),1).*sqrt(-1).*imag(z);



function[]=ellipseplot_fig

a=ones(10,1);
b=(1:10)'./10;
[k,l]=ab2kl(a,b);
th=linspace(0,pi,10)';
x=linspace(0,1,10)'*20;
ellipseplot(k,l,th,x,'npoints',16,'g')
title('Counterclockwise rotating ellipse becoming circle')
set(gca,'dataaspectratio',[1 1 1])
