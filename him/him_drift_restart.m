function [ array ] = him_drift_restart( inpath )

files=dir([inpath filesep 'HIM.res_*.nc']);
fill4d=zeros(length(files),49,210,360);

array.ptemp=fill4d;
array.sal=fill4d;
array.h=fill4d;
array.uh=fill4d;
array.vh=fill4d;
tic;
for i=1:length(files)
    disp(sprintf('File %d/%d Elapsed Time: %f',i,length(files),toc));
    ncfile=[inpath filesep files(i).name];
    array.ptemp(i,:,:,:)=nc_varget(ncfile,'Temp');
    array.sal(i,:,:,:)=nc_varget(ncfile,'Salt');
    array.h(i,:,:,:)=nc_varget(ncfile,'h');
    array.uh(i,:,:,:)=nc_varget(ncfile,'uh');
    array.vh(i,:,:,:)=nc_varget(ncfile,'vh');
    
end