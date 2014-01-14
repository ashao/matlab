function [ ] = autoguess_all_aao( )

lags=-60:5:60;
lags(lags==0)=[];

for ii=1:length(lags)    
    guesspath=sprintf('/scratch/data/jetmodel/toprocess/aao/lag_00/');
    inpath=sprintf('/scratch/data/jetmodel/aao/lag_%02d/',lags(ii));
    outpath=sprintf('/scratch/data/jetmodel/toprocess/aao/lag_%02d/',lags(ii));
    auto_lag_endpoints(guesspath,inpath,outpath)
    
end