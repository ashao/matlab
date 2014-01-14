function [ atmhist ] = atmospheric_histories( )

    % Download most recent 
    sourcelink='http://cdiac.ornl.gov/ftp/oceans/CFC_ATM_Hist/CFC_Atm_Hist_2011.csv';
    filename=sprintf('CFC_Atm_Hist_2011.%d.csv',floor(now));
    urlwrite(sourcelink,filename);
    data=csvread(filename,1,0);
    
    atmhist.time=data(:,1);
    atmhist.cfc11.nval=data(:,2);
    atmhist.cfc11.sval=data(:,3);
    
    atmhist.cfc12.nval=data(:,4);
    atmhist.cfc12.sval=data(:,5);

    atmhist.sf6.nval=data(:,10);
    atmhist.sf6.sval=data(:,11);
    delete(filename)
end