function [gh11,nh11,gh12,nh12,gh22,nh22,gh13,nh13,gh23,nh23,gh33,nh33]=variof2(x1,x2,x3,icode);
%
% function [gh11,nh11,gh12,nh12,gh22,nh22,gh13,nh13,gh23,nh23,gh33,nh33]=variof2(x1,x2,x3,icode);
%
% function to compute variograms, cross-variograms, covariograms,
% cross-covariograms and pseudo-cross-variograms in 1D or 2D for up to 3 variables.
% the data are on a (possibly incomplete) regular grid
% the program computes variograms in the frequency domain using 2D-FFT.
%
% input: x1,x2,x3:   data matrices. Each matrix is either empty or of size n x m
%                    Missing values are indicated by NaN
%        icode: a code to indicate which function to compute
%               =1 : variograms and cross-variograms;
%               =2 : covariograms and cross-covariograms
%               =3 : variograms and pseudo-cross-variograms
%               =4 : on calcule une moyenne pour l'ensemble du champ au lieu de selon les lags
%                    (transl: We calculate an average for the whole field instead of for the lags, M. Hood)
%
% output: ghij: (direct- and cross-) variograms for variables i and j, or covariograms or
%               pseudo-variograms depending on icode.
%         nhij: number of pairs available to compute the structural function
%
% this program uses the functions FFT2, IFFT2, FFT2SHIFT and CONJ which are
% standard MATLAB functions.
% author: D. Marcotte, dmarcotte@mail.polymtl.ca
gh11=[];nh11=[];gh12=[];nh12=[];gh22=[];nh22=[];gh13=[];nh13=[];gh23=[];nh23=[];gh33=[];nh33=[];
[n,p]=size(x1);				% dimensions of data matrix
[n2,p2]=size(x2);
[n3,p3]=size(x3);

% find the closest multiple of 8 to obtain a good compromise between
% speed (a power of 2) and memory required

nrows=2*n-1;
ncols=2*p-1;
nr2=ceil(nrows/8)*8;
nc2=ceil(ncols/8)*8;

% form an indicator  matrix: 1's for all data values
%                            0's for missing values
% in data matrix, replace missing values by 0;

x1id=~isnan(x1);			% 1 for a data value; 0 for missing
x1(~x1id)=zeros(sum(sum(~x1id)),1);	% missing replaced by 0

% compute fourier transforms
if icode==4;
  m1=sum(sum(x1(x1id)))/sum(sum(x1id));
  x1(x1id)=x1(x1id)-m1;
end

fx1=fft2(x1,nr2,nc2);			% fourier transform ox x1
fx1_x1=fft2(x1.*x1,nr2,nc2);		% fourier transform of x1*x1
fx1id=fft2(x1id,nr2,nc2);
nh11=round(real(ifft2(conj(fx1id).*fx1id))); % number of pairs for x1 variogram

% do the same with x2 (if defined)
% compute number of pairs for cross-structural functions selected (icode)

if n2>0;				% begin test on n2
  x2id=~isnan(x2);
  x2(~x2id)=zeros(sum(sum(~x2id)),1);
  if icode==4;
    m2=sum(sum(x2(x2id)))/sum(sum(x2id));
    x2(x2id)=x2(x2id)-m2;
  end

  fx2=fft2(x2,nr2,nc2);
  fx2_x2=fft2(x2.*x2,nr2,nc2);
  fx2id=fft2(x2id,nr2,nc2);
  nh22=round(real(ifft2(conj(fx2id).*fx2id)));

  % number of pairs for cross-structural function x1-x2

  if icode==1
    fx12id=fft2(x1id.*x2id,nr2,nc2);
    nh12=round(real(ifft2(conj(fx12id).*fx12id)));
  else
    nh12=round(real(ifft2(conj(fx1id).*fx2id)));
  end
end					% end test on n2

% do the same with x3 (if defined)
% compute number of pairs for cross-structural functions selected (icode)

if n3>0;				% begin test on n3
  x3id=~isnan(x3);
  x3(~x3id)=zeros(sum(sum(~x3id)),1);
  if icode==4;
    m3=sum(sum(x3(x3id)))/sum(sum(x3id));
    x3(x3id)=x3(x3id)-m3;
  end

  fx3=fft2(x3,nr2,nc2);
  fx3_x3=fft2(x3.*x3,nr2,nc2);
  fx3id=fft2(x3id,nr2,nc2);
  nh33=round(real(ifft2(conj(fx3id).*fx3id)));

  % number of pairs for cross-structural function x1-x3, x2-x3

  if icode==1
    fx13id=fft2(x1id.*x3id,nr2,nc2);
    nh13=round(real(ifft2(conj(fx13id).*fx13id)));
    fx23id=fft2(x2id.*x3id,nr2,nc2);
    nh23=round(real(ifft2(conj(fx23id).*fx23id)));
  else
    clear x1 x2 x3 x1id x2id x3id
    nh13=round(real(ifft2(conj(fx1id).*fx3id)));
    nh23=round(real(ifft2(conj(fx2id).*fx3id)));
  end

end					% end test on n3

% compute the different structural functions according to icode

if icode==1;				% variograms and cross-variograms are computed

  if n3>0;
    gh33=real(ifft2(conj(fx3id).*fx3_x3+conj(fx3_x3).*fx3id-2*conj(fx3).*fx3))./max(nh33,1)/2;
    clear fx3id fx3_x3 fx3

    t1=fft2(x1.*x3id,nr2,nc2);
    t2=fft2(x3.*x1id,nr2,nc2);
    t12=fft2(x1.*x3,nr2,nc2);
    gh13=real(ifft2(conj(fx13id).*t12+conj(t12).*fx13id-conj(t1).*t2-t1.*conj(t2)))./max(nh13,1)/2;
    clear fx13id

    t1=fft2(x2.*x3id,nr2,nc2);
    t2=fft2(x3.*x2id,nr2,nc2);
    t12=fft2(x2.*x3,nr2,nc2);
    clear x3 x3id
    gh23=real(ifft2(conj(fx23id).*t12+conj(t12).*fx23id-conj(t1).*t2-t1.*conj(t2)))./max(nh23,1)/2;
    clear fx23id
  end

  if n2>0;
    gh22=real(ifft2(conj(fx2id).*fx2_x2+conj(fx2_x2).*fx2id-2*conj(fx2).*fx2))./max(nh22,1)/2;
    clear fx2id fx2_x2 fx2

    t1=fft2(x1.*x2id,nr2,nc2);		% FFT on points available for cross-variogram computation
    clear x2id
    t2=fft2(x2.*x1id,nr2,nc2);
    clear x1id
    t12=fft2(x1.*x2,nr2,nc2);
    clear x1 x2
    gh12=real(ifft2(conj(fx12id).*t12+conj(t12).*(fx12id)-conj(t1).*t2-t1.*conj(t2)))./max(nh12,1)/2;
    clear fx12id t1 t2 t12
  end

  gh11=real(ifft2(conj(fx1id).*fx1_x1+conj(fx1_x1).*fx1id-2*conj(fx1).*fx1))./max(nh11,1)/2;

elseif icode==2;			% covariograms and cross-covariograms are computed

  clear fx1_x1 fx2_x2 fx3_x3

  if n3>0;
    m1=real(ifft2(conj(fx3id).*fx3))./max(nh33,1); % computes the tail x3 means
    m2=real(ifft2(fx3id.*conj(fx3)))./max(nh33,1); % computes the head x3 means
    gh33=real((ifft2(conj(fx3).*fx3))./max(nh33,1)-m1.*m2);

    m1=real(ifft2(conj(fx1).*fx3id))./max(nh13,1); % computes the tail x1 means
    m2=real(ifft2(conj(fx1id).*fx3))./max(nh13,1); % computes the head x3 means
    gh13=real((ifft2(conj(fx1).*fx3))./max(nh13,1)-m1.*m2);

    m1=real(ifft2(conj(fx2).*fx3id))./max(nh23,1); % computes the tail x2 means
    m2=real(ifft2(conj(fx2id).*fx3))./max(nh23,1); % computes the head x3 means
    gh23=real((ifft2(conj(fx2).*fx3))./max(nh23,1)-m1.*m2);
    clear fx3 fx3id
  end

  if n2>0;
    m1=real(ifft2(conj(fx2).*fx2id))./max(nh22,1); % computes the tail x2 means
    m2=real(ifft2(conj(fx2id).*fx2))./max(nh22,1); % computes the head x2 means
    gh22=real((ifft2(conj(fx2).*fx2))./max(nh22,1)-m1.*m2);

    m1=real(ifft2(conj(fx1).*fx2id))./max(nh12,1); % computes the tail x1 means
    m2=real(ifft2(conj(fx1id).*fx2))./max(nh12,1); % computes the head x2 means
    gh12=real(ifft2(conj(fx1).*fx2))./max(nh12,1)-m1.*m2;
    clear fx2 fx2id
  end

  m1=real(ifft2(conj(fx1).*fx1id))./max(nh11,1);   % computes the tail x1 means
  m2=real(ifft2(conj(fx1id).*fx1))./max(nh11,1);   % computes the head x1 means
  gh11=real((ifft2(conj(fx1).*fx1))./max(nh11,1)-m1.*m2);

elseif icode==3				% variograms and pseudo-cross-variograms are computed

  if n3>0;
    gh33=real(ifft2(conj(fx3id).*fx3_x3+conj(fx3_x3).*fx3id-2*conj(fx3).*fx3))./max(nh33,1)/2;
    gh13=real(ifft2(fx3id.*conj(fx1_x1)+conj(fx1id).*fx3_x3-2*conj(fx1).*fx3))./max(nh13,1)/2;
    gh23=real(ifft2(fx3id.*conj(fx2_x2)+conj(fx2id).*fx3_x3-2*conj(fx2).*fx3))./max(nh23,1)/2;
    clear fx3id fx3 fx3_x3
  end

  if n2>0;
    gh22=real(ifft2(conj(fx2id).*fx2_x2+conj(fx2_x2).*fx2id-2*conj(fx2).*fx2))./max(nh22,1)/2;
    gh12=real(ifft2(fx2id.*conj(fx1_x1)+conj(fx1id).*fx2_x2-2*conj(fx1).*fx2))./max(nh12,1)/2;
    clear fx2id fx2 fx2_x2
  end

  gh11=real(ifft2(conj(fx1id).*fx1_x1+conj(fx1_x1).*fx1id-2*conj(fx1).*fx1))./max(nh11,1)/2;

elseif icode==4;
  clear fx1_x1 fx2_x2 fx3_x3

  if n3>0;
    gh33=real((ifft2(conj(fx3).*fx3))./max(nh33,1));
    gh13=real((ifft2(conj(fx1).*fx3))./max(nh13,1));
    gh23=real((ifft2(conj(fx2).*fx3))./max(nh23,1));
    clear fx3 fx3id
  end

  if n2>0;
    gh22=real(ifft2(conj(fx2).*fx2))./max(nh22,1);
    gh12=real(ifft2(conj(fx1).*fx2))./max(nh12,1);
    clear fx2 fx2id
  end

  gh11=real((ifft2(conj(fx1).*fx1))./max(nh11,1));

end					% end test on icode

% reduce matrices to required size

nh11=[nh11(1:n,1:p) nh11(1:n,nc2-p+2:nc2);nh11(nr2-n+2:nr2,1:p) nh11(nr2-n+2:nr2,nc2-p+2:nc2)];
gh11=[gh11(1:n,1:p) gh11(1:n,nc2-p+2:nc2);gh11(nr2-n+2:nr2,1:p) gh11(nr2-n+2:nr2,nc2-p+2:nc2)];

if n2>0
  nh22=[nh22(1:n,1:p) nh22(1:n,nc2-p+2:nc2);nh22(nr2-n+2:nr2,1:p) nh22(nr2-n+2:nr2,nc2-p+2:nc2)];
  gh22=[gh22(1:n,1:p) gh22(1:n,nc2-p+2:nc2);gh22(nr2-n+2:nr2,1:p) gh22(nr2-n+2:nr2,nc2-p+2:nc2)];
  nh12=[nh12(1:n,1:p) nh12(1:n,nc2-p+2:nc2);nh12(nr2-n+2:nr2,1:p) nh12(nr2-n+2:nr2,nc2-p+2:nc2)];
  gh12=[gh12(1:n,1:p) gh12(1:n,nc2-p+2:nc2);gh12(nr2-n+2:nr2,1:p) gh12(nr2-n+2:nr2,nc2-p+2:nc2)];

  if n3>0
    nh33=[nh33(1:n,1:p) nh33(1:n,nc2-p+2:nc2);nh33(nr2-n+2:nr2,1:p) nh33(nr2-n+2:nr2,nc2-p+2:nc2)];
    gh33=[gh33(1:n,1:p) gh33(1:n,nc2-p+2:nc2);gh33(nr2-n+2:nr2,1:p) gh33(nr2-n+2:nr2,nc2-p+2:nc2)];
    nh13=[nh13(1:n,1:p) nh13(1:n,nc2-p+2:nc2);nh13(nr2-n+2:nr2,1:p) nh13(nr2-n+2:nr2,nc2-p+2:nc2)];
    gh13=[gh13(1:n,1:p) gh13(1:n,nc2-p+2:nc2);gh13(nr2-n+2:nr2,1:p) gh13(nr2-n+2:nr2,nc2-p+2:nc2)];
    nh23=[nh23(1:n,1:p) nh23(1:n,nc2-p+2:nc2);nh23(nr2-n+2:nr2,1:p) nh23(nr2-n+2:nr2,nc2-p+2:nc2)];
    gh23=[gh23(1:n,1:p) gh23(1:n,nc2-p+2:nc2);gh23(nr2-n+2:nr2,1:p) gh23(nr2-n+2:nr2,nc2-p+2:nc2)];
  end
end

% shift all the matrices so that the 0 lag appears at the center of each matrix

nh11=fftshift(nh11);
gh11=fftshift(gh11);
nh22=fftshift(nh22);
gh22=fftshift(gh22);
nh12=fftshift(nh12);
gh12=fftshift(gh12);
nh33=fftshift(nh33);
gh33=fftshift(gh33);
nh13=fftshift(nh13);
gh13=fftshift(gh13);
nh23=fftshift(nh23);
gh23=fftshift(gh23);
