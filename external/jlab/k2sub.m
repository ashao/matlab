function[m,n,index]=k2sub(K,dk,N)
%K2SUB Convert a complex-valued wavenumber into an I,J subscript pair
%  
%   [I,J,INDEX]=K2SUB(K,DK,N), K is a complex-valued wavenumber which
%   occurs on an NxN wavenumber grid with wavenumber spacing DK,
%   returns the I,J indices corresponding to the location of K on the
%   grid.  
%  
%   Only I,J values >0 and <N+1 are returned, and these correspond to   
%   locations INDEX in K(:).  
%  
%   K2SUB is inverted by SUB2K.  See also WAVEGRID.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003, 2004 J.M. Lilly --- type 'help jlab_license' for details        

if strcmp(K, '--t')
  k2sub_test,return
end

n=real(K)./dk+(N+1)/2;
m=imag(K)./dk+(N+1)/2;

n=round(n);
m=round(m);

index=find(n>0 & n<N+1 & m>0 & m<N+1);
[n,m]=vindex(n(:),m(:),index,1);


function[]=k2sub_test

k=sub2k((1:5),0*(1:5)+2+1,1,5);
[m,n]=k2sub(k,1,5);
bool=aresame(m,(1:5)') && aresame(n,3+0*n);
reporttest('K2SUB for 1 x 5 unit spacing wavegrid',bool)

k=sub2k((1:4),0*(1:4)+2,1,4);
[m,n]=k2sub(k,1,4);
bool=aresame(m,(1:4)') && aresame(n,2+0*n);
reporttest('K2SUB for 1 x 5 unit spacing wavegrid',bool)

dk=1;N=5;
[ii,jj]=k2sub(wavegrid(dk,N),dk,N);
[ii2,jj2]=ind2sub([5 5],(1:25));
bool=allall(ii(:)==ii2(:)) && allall(jj(:)==jj2(:));
reporttest('K2SUB for 5 x 5 unit spacing wavegrid',bool)
  
