function[ir,jr,wr,fr,br,cr]=ridgewalk(varargin)
% RIDGEWALK  Extract wavelet transform ridges, including bias estimates.
%  
%   [IR,JR,XR,FR]=RIDGEWALK(W,FS) where W is a wavelet transform matrix
%   at frequecies FS, such as that returned by MORSEWAVE, returns the 
%   wavelet ridges of transform W.
%
%   The columns of W correspond to different frequencies, specified by the
%   frequency array FS, at which the wavelet transform was performed.  Note
%   that FS assumes a unit sample rate.  
%
%   The frequencies FS are expected to be ordered from highest to lowest.  
%
%   RIDGEWALK returns the following quantities along ridges
%
%       IR     Ridge indices into rows of W (time) 
%       JR     Ridge indices into columns of W (scale) 
%       XR     Estimated analtic signal along the ridge
%       FR     Transform frequency values in radian frequency
%
%   All output variables are vectors of the same length.  A NAN appears in
%   each value at the bottom of every ridge to separate different ridges. 
%
%   RIDGEWALK(DT,W,FS,...) uses a sample rate DT to compute the ridge
%   frequency FR.  The default value of DT is unity.  This does not affect
%   the specification of FS, which is given in terms of a unit sample rate.
%   _______________________________________________________________________
%
%   Options
%
%   RIDGEWALK(...,{N,CHI,ALG}) specifies options for the ridge computation.
%
%        N  -- Removes all ridges of less than N periods in length
%      CHI  -- Removes all small amplitude ridge points having |W|<CHI
%      ALG  -- String specifying algorithm; see below
%
%   RIDGEWALK(...,{ 1.5,   0   ,'amp'}) is the default.
%                     |    |      | 
%                     N   CHI    ALG
%
%   See below for details on the ridge algorithm type.
%   _______________________________________________________________________
%
%   Ridge types
%
%   Two different definitions may be used to locate the ridges.
%
%   ALG determines the ridge type, where ALG may be either:
%
%         'phase'       Rate of transform change of phase definition
%         'amplitude'   Maxima of transfom amplitude definition
% 
%   If ALG is not specified, 'amplitude' is used by default.
%
%   In practice, these usually do not differ much from one another.  An
%   examination of the difference between phase and amplitude ridges may be
%   found in 
%
%      Lilly and Olhede (2010).  On the analytic wavelet transform.
%
%   For reasons given therein, we prefer amplitude ridges.
%   _______________________________________________________________________
%
%   Time-dependent frequency range
%
%   The ridge curves may be limited to a time-varying frequency range.
%
%   RIDGEWALK(...,{FMAX,FMIN,N,CHI,ALG}) additionally specifies a 
%   maximim frequency FMAX and minumum frequency FMIN for the ridges.  
%   Only ridge points between these two frequencies are returned.
%
%   FMAX and FMIN are both *radian* frequencies per unit time as specified
%   by DT. DT is optional and its default value is unity.  Thus FMAX and 
%   FMIN are directly comparable to the ridge frequency FR. 
%   
%   Both FMAX and FMIN are the same length as W. 
%   _______________________________________________________________________
%
%   Optional output  
%
%   [IR,JR,XR,FR,BR,CR]=RIDGEWALK(...) optionally outputs two additional
%   quantities along the ridges.
%
%       BR     Instantaneous bandwidth  
%       CR     Instantaneous curvature  
%
%   When the 'bias parameters' BR and CR are small, the signal is 
%   accurately estimated, as discussed in 
%
%      Lilly and Olhede (2010).  On the analytic wavelet transform.
%
%   Note that in Lilly and Olhede (2010), the relative bandwidth BR is
%   denoted upsilon/omega, while the complex-valued chirp rate CR is 
%   d/dt [omega-i*upsilon] / omega^2.
%   _______________________________________________________________________
%
%   Joint ridges
%
%   RIDGEWALK(W1,W2,...,WM,FS) finds the joint ridges of M transforms.
%   These transforms all have the same size.  
%
%   In this case, there is only one set of ridges but M different values. 
%   Thus all output arguments except IR and JR are two-dimensional arrays
%   with M elements in the second dimension.
%
%   RIDGEWALK(W,FS) with W a 3-D array having M ``pages'' such that
%   W(:,:,1)=W1, W(:,:,2)=W2,...,W(:,:,M)=WM, also works.
%
%   For details on joint ridges, see
%
%      Lilly and Olhede (2011).  Theory and Analyisis of Modulated 
%            Multivariate Oscillations, Part I: Fundamentals
%
%   For joint ridges, FR, BR, and CR remain the instantaneous moments of 
%   the individual signals, and thus have the same size as XR.  To compute
%   the joint instantaneous moments from these, use JOINTFREQ.
%   _______________________________________________________________________
%
%
%   Interscale interpolation
%   
%   RIDGEWALK interpolates among discrete scale levels to yield more
%   accurate values for the ridge quantities XR and FR using a fast
%   quadratic interpolation.  
%   
%   See RIDGEINTERP and QUADINTERP for details.
%   _______________________________________________________________________
%
%   See also WAVETRANS, RIDGEMAP, RIDGEINTERP.
%
%   'ridgewalk --t' runs a test.
%   'ridgewalk --f' generates a sample figure.
%
%   Usage: [ir,jr,wr,fr]=ridgewalk(w,fs);
%          [ir,jr,wr,fr]=ridgewalk(w,fs,{N,CHI,'amp'});
%          [ir,jr,wr,fr]=ridgewalk(dt,w,fs,{N,CHI,'amp'});
%          [ir,jr,wr,fr]=ridgewalk(dt,w,fs,{1.5,0,'amp'});
%   _______________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details
%   


%   RIDGEWALK can be used to create a de-biased estimate of the signal,
%   following Lilly and Olhede (2010).  This estimate is given by
%
%       XR_DB= XR - (XR/2)*(BR^2+i*CR)*P^2
%
%   where XR, BR, and CR are defined above, where P is the wavelet 
%   time-frequency product, and i=SQRT(-1).  For the generalized Morse 
%   wavelets, P=SQRT(BETA*GAMMA).



%   ALPHA  -- Controls agressiveness of chaining across scales.
%
%   Chaining parameter
%
%   The chaining parameter ALPHA specifies the agressiveness with which
%   ridge points are chained across scales.  
%
%   The default value of ALPHA is one-half.
%
%   Increase ALPHA to chain ridges more agressively across scales, or 
%   decrease ALPHA to supress chaining across scales.
%   ___________________________________________________________________
%
%   ALPHA is defined as a normalized frequency difference
%
%         ALPHA  =  DOMEGA / OMEGA
%
%   where OMEGA is the transform frequency, and DOMEGA is the difference
%   between the frequency predicted for the next point based on the
%   transform at a "tail", and the actual frequency at prospective "heads".  
%
%   The chaining parameter is defined in such a way that it does not
%   need to be changed as time sampling or frequency sampling changes.
%   However, for strongly chirping signals or weakly chirping, noisy 
%   signals, better performance may be obtianed by adjusting it.
%   ______________________________________________________________
%
%   Spurious ridge points
%
%   RIDGEWALK has a continigency for rejecting spurious ridge points.
%   These tend to occur on the flanks of interesting signals, and 
%   reflect the wavelet structure rather than the signal structure.
%
%   See ISRIDGEPOINT for details.


if strcmp(varargin{1}, '--t')
    ridgewalk_test,return
elseif strcmp(varargin{1}, '--f')
    ridgewalk_figure,return
end

if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else
    dt=1;
end

if iscell(varargin{end})
    params=varargin{end};
    varargin=varargin(1:end-1);
else
    params={1.5,0,1/4,'amp'};
end

fs=varargin{end};
if ~aresame(sort(fs(:)),flipud(fs(:)))
    error('The frequencies FS should be sorted from highest to lowest.')
end
varargin=varargin(1:end-1);


%/********************************************************************
%Sorting out input params
if length(params{1})>1
    fmax=row2col(params{1}).*dt;
    fmin=row2col(params{2}).*dt;
    params=params(3:end);
else
    fmax=[];
    fmin=[];
end
    
if length(params)==3
    alpha=1/4;
    N=params{1};
    chi=params{2};
    alg=params{3};
else
    N=params{1};
    chi=params{2};
    alpha=params{3};
    alg=params{4};
end
%\********************************************************************

if length(varargin)>1
    w=zeros([size(varargin{1},1) size(varargin{1},2) length(varargin)]);
    for i=1:length(varargin);
        w(:,:,i)=varargin{i};
    end
else
    w=varargin{1};
end
if size(w,3)
    disp(['RIDGEWALK detecting a set of ' int2str(size(w,3)) ' transforms.'])
end

%Changes below to output deviation vectors
if nargout<=4
    om=instfreq(dt,w,'endpoint');
elseif nargout==5
    [om,up]=instfreq(dt,w,'endpoint');
elseif nargout==6
    [om,up,curv]=instfreq(dt,w,'endpoint');
end

[bool,rq,wjoint,omjoint]=isridgepoint(w,fs,chi,alg,fmin,fmax);  

disp('RIDGEWALK chaining ridges.')
[id,ir,jr,wr,fr]=ridgechains(fs,N,bool,wjoint,omjoint,alpha);
if ~isempty(id)
     [id,ir,jr,wr,fr]=colbreaks(id,ir,jr,wr,fr);
     [wr,fr]=ridgeinterp(wjoint,fs,rq,ir,jr,w,om);
end

if nargout ==5
    br=ridgeinterp(wjoint,fs,rq,ir,jr,up);
elseif nargout ==6
    [br,cr]=ridgeinterp(wjoint,fs,rq,ir,jr,up,curv);
end

disp('RIDGEWALK finished.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[bool,rq,w,om]=isridgepoint(w,fs,chi,str,fmin,fmax)
%ISRIDGEPOINT  Finds wavelet ridge points using one of several criterion.
%
%   BOOL=ISRIDGEPOINT(W,FS,CHI,STR) where W is a wavelet transform matrix at 
%   *cyclic* frequecies FS, finds all ridge points of W with amplitudes
%   |W| exceeding the amplitude cutoff A.  Several different different 
%   ridge defintions may be used and are specified by STR.
%
%   BOOL is a matrix of the same size as W, which is equal to one for 
%   those elements of W which are ridge points, and zero otherwise.
%
%   STR may be either of the following:
%
%        'phase'       Rate of transform change of phase definition
%        'amplitude'   Maxima of transfom amplitude definition
%
%   For all definitions, ISRIDGEPOINT rejects spurious ridge points.
%   These tend to occur on the flanks of interesting signals, and 
%   reflect the wavelet structure rather than the signal structure.
%
%   A ridge point is considered spurious if either it is located at an
%   amplitude minima, or if the frequency anomaly (transform frequency
%   minus scale frequency) is a maximum.
%
%   See also RIDGEQUANTITY, RIDGEWALK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2009 J.M. Lilly --- type 'help jlab_license' for details
 
%        'groove'      Joint amplitude / phase definition

disp('RIDGEWALK looking for ridge points...')

[rq,om]=ridgequantity(w,fs,str);
if size(w,3)~=1
    phaseavg=frac(sum(abs(w).*w,3),sum(abs(w).^2,3));
    w=sqrt(sum(abs(w).^2,3)).*rot(angle(phaseavg));
end

rqm=vshift(rq,-1,2);
rqp=vshift(rq,+1,2);

%This is d/ds < 0 since scale decreases in columns
bool=(rqm<0&rqp>=0)|(rqm<=0&rqp>0);
   
err=abs(rq);

%Ensure maximum not minimum
bool((bool&vshift(bool,1,2))&err>vshift(err,1,2))=0; 
bool((bool&vshift(bool,-1,2))&err>vshift(err,-1,2))=0; 

bool1= ~isnan(w);   %Remove NANs
bool2=~(abs(w)<chi);  %Remove those less than cutoff amplitude
bool=bool.*bool1.*bool2;
bool(:,[1 end])=0;

%Running FMIN and FMAX 
if ~isempty(fmin)
    fsmat=vrep(col2row(fs(:)),size(w,1),1);
    fmin=vrep(fmin,size(w,2),2);
    fmax=vrep(fmax,size(w,2),2);
    bool3=(fsmat>=fmin)&(fsmat<=fmax);
    bool=bool.*bool3;
end

disp(['RIDGEWALK found ' int2str(length(find(bool))) ' ridge points.'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[id,ii,jj,xr,fr]=ridgechains(fs,N,bool,x,f,alpha)
%RIDGECHAINS  Forms ridge curves by connecting transform ridge points.
%
%   [ID,IR,JR,XR]=RIDGECHAINS(N,BOOL,W) forms chains of ridge points
%   of wavelet transform W.  
%
%   Ridge points are all points of W for which BOOL, a matrix of the 
%   same size as W, equals one.  Only ridges of at least N periods in 
%   length are returned.
%   
%   ID is a unique ID number assigned to each ridge.  IR and JR are 
%   the time- and scale-indices along the ridges.  XR is the wavelet
%   transform along the ridge. 
% 
%   All output variables are the same size.
% 
%   See also RIDGEWALK.
%
%   Usage:  [id,ii,jj,xr]=ridgechains(N,bool,x);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2007 J.M. Lilly --- type 'help jlab_license' for details
 
if isempty(find(bool,1))
    [id,ii,jj,xr,fr]=vempty;return
end
    
dfdt=vdiff(f,1);

indexridge=find(bool);
[ii,jj]=ind2sub(size(bool),indexridge);
[ii,sorter]=sort(ii);
vindex(jj,indexridge,sorter,1);

%Using new algorithm as of November 2007, faster and also prevents ridge breaking

xr=x(indexridge);             %Transform value along ridge
fr=f(indexridge);             %Frequency along ridge
fsr=fs(jj);                   %Scale frequency along ridge
fr_next=fr+dfdt(indexridge);  %Predicted frequency at next point
fr_prev=fr-dfdt(indexridge);  %Predicted frequency at previous point

%figure,plot(fr,'.')
%figure,plot(ii,jj,'r.'),figure,plot(ii,abs(xr),'r.')
cumbool=cumsum(bool,2);
J=maxmax(cumbool(:,end));

[indexmat,nextindexmat,iimat,jjmat,fsmat,frmat,fr_nextmat,fr_prevmat]=vzeros(size(f,1),J,'nan');

%Indices for this point
indexmat(sub2ind(size(indexmat),ii,cumbool(indexridge)))=(1:length(ii));
nonnanindex=find(~isnan(indexmat));

%Don't overwrite the original variables
[ii1,jj1,fsr1,fr1,fr_next1,fr_prev1]=vindex(ii,jj,fsr,fr,fr_next,fr_prev,indexmat(nonnanindex),1);
vindexinto(iimat,jjmat,fsmat,frmat,fr_nextmat,fr_prevmat,ii1,jj1,fsr1,fr1,fr_next1,fr_prev1,nonnanindex,0);
clear ii1 jj1 fsr1 fr1 fr_next1 fr_prev1

%Time difference from points here to next points
dii=vsum(vshift(iimat,1,1)-iimat,2);  %To get rid of nans

clear iimat

%Scale frequency difference from points here to next points
fsmat3=vrep(fsmat,J,3);
frmat3=vrep(frmat,J,3);

%Predicted minus actual frequency at this point
fr_nextmat3=vrep(fr_nextmat,J,3);
df1=permute(vshift(frmat3,1,1),[1 3 2])-fr_nextmat3;
df1=frac(df1,frmat3);

clear fr_nextmat3

%Expected minus actual frequency at next point
fr_prevmat3=vrep(fr_prevmat,J,3);
df2=permute(vshift(fr_prevmat3,1,1),[1 3 2])-frmat3;
df2=frac(df2,frmat3);

df=frac(abs(df1)+abs(df2),2);
%size(ii)
%figure,plot(df(:),'.'),ylog

df(df>alpha)=nan;

clear fr_prevmat3 fr_mat3 df1 df2 


%Keep when they are the same 
%Set df to nan except when min along one direction
[mindf,jjmin]=min(df,[],2);
iimin=vrep((1:size(df,1))',size(df,2),2);
kkmin=vrep(1:size(df,2),size(df,1),1);
df=nan*df;
df(sub2ind(size(df),iimin,squeeze(jjmin),kkmin))=mindf;
clear iimin jjmin kkmin 

[mindf,jjmin]=min(df,[],3);
index=find(~isnan(mindf));

if ~isempty(index)
    mindf=mindf(index);
    [ii2,jj2]=ind2sub(size(indexmat),index);
    [ii2,jj2,index,mindf]=vindex(ii2,jj2,index,mindf,find(ii2<size(x,1)),1);
    index2=sub2ind(size(indexmat),ii2+1,jjmin(index));
    nextindexmat(index)=indexmat(index2);
end

%df=nan*ii;df(ii2)=mindf;

id=nan*ii;
%Assign a unique number to all ridge points
for i=1:size(nextindexmat,1)
    id(nonnan(indexmat(i,:)))=nonnan(indexmat(i,:));
end

%Reassign number for linked points
for i=1:size(nextindexmat,1)
    id(nonnan(nextindexmat(i,:)))=id(indexmat(i,~isnan(nextindexmat(i,:))));
end
[id,sorter]=sort(id);
vindex(ii,jj,indexridge,xr,fr,sorter,1);

%Remove ridge lines of length shorter than a specified length
lr=ridgelen(id,ii,jj,fr);
vindex(id,ii,jj,indexridge,xr,fr,find(lr>=N),1); 
       
disp(['RIDGEWALK pruning to ' int2str(length(id)) ' ridge points.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Tests and figures below here


function[]=ridgewalk_figure
 
%/******************************************************************
load ebasn
use ebasn


len=cellength(lat);
index=find(len>200);
lato=30;

id=id(index);num=num(index);lat=lat(index);lon=lon(index);
p=p(index);t=t(index);s=s(index);

index=24;
id=id{index};num=num{index};lat=lat{index};lon=lon{index};
p=p{index};t=t{index};s=s{index};
dt=num(2)-num(1);

ga=3;be=3;
fs=morsespace(ga,be,{0.05,2*pi/3},2*pi/100,8);


%Compute wavelet transforms using generalized Morse wavelets

cx=fillbad(latlon2xy(lat,lon,30,-25));
cv=latlon2uv(num,lat,lon);

wx=wavetrans(real(cx),{1,ga,be,fs,'bandpass'},'mirror');
wy=wavetrans(imag(cx),{1,ga,be,fs,'bandpass'},'mirror');

[ir,jr,xr,fr]=ridgewalk(dt,wx,wy,fs,{2*morseprops(ga,be),0,'amp'});  
[xhat,fhat]=ridgemap([length(cx) 2],xr,fr,ir);
fbar=instfreq(xhat,2);
    

ci=(0:5:65);
numo=datenum(1986,1,1)-1;
[h,hl]=wavespecplot(num-numo,cv,2*pi./fs,sqrt(abs(wx).^2+abs(wy).^2),1,ci);
linestyle -h hl k k--
axes(h(1)),ylim([-18 18]),ylabel('Current Speed (cm/s)'),title('Bivariate Ridge Method Example')
text(-90,15,'(a)')

axes(h(2)),caxis([0 40]),colormap gray,flipmap,ylim([3.6 60]),hold on
plot(num-numo,2*pi./fbar,'w','linewidth',4)
plot(num-numo,2*pi./fbar,'k','linewidth',2)

xlabel('Day of Year 1986'),ylabel('Period in Days')
set(gca,'ytick',2.^(2:.5:5.5))
set(gca,'yticklabel',[' 4';'  ';' 8';'  ';'16';'  ';'32';'  '])
inticks
text(-90,4.5,'(b)')

orient landscape
fontsize 12 10 10 10
set(gcf,'paperposition',[1 1 9 5])


function[]=ridgewalk_test

load npg2006
use npg2006

%Decide on frequencies
fs=2*pi./(logspace(log10(10),log10(100),50)');

%Compute wavelet transforms using generalized Morse wavelets
[wx,wy]=wavetrans(real(cx),imag(cx),{1,2,4,fs,'bandpass'},'mirror');
[wp,wn]=vectmult(tmat,wx,wy);


%Form ridges of component time series
[ir,jr,wr,fr]=ridgewalk(dt,wn,fs,{1.5,0,'phase'}); 
[ir2,jr2,wr2,fr2]=ridgewalk(dt,wn,fs,{1.5,0,'amplitude'});

err=vsum(abs(wr-wr2).^2,1)./vsum(abs(wr).^2,1);
reporttest('RIDGEWALK phase and amplitude signal estimation error test for NPG-06',err<1e-3)

[ir,jr,wr,fr]=ridgewalk(dt,wp,wn,fs,{3,0,'amp'});   
[ir2,jr2,wr2,fr2]=ridgewalk(dt,wx,wy,fs,{3,0,'amp'});   
[wrp2,wrn2]=vectmult(tmat,wr2(:,1),wr2(:,2));


err=vmean((abs(wr-[wrp2 wrn2])./sqrt(abs(wr).^2+abs([wrp2 wrn2])).^2),1);
reporttest('RIDGEWALK XY vs. PN invariance joint ridge test for NPG-06',allall(err<1e-10))







