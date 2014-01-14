function [ ] = him_proc_benchmark( ncfile, field, nprocs, months_per_proc )

if matlabpool('size')==0
    matlabpool(nprocs);
end

if matlabpool('size')~=nprocs
    matlabpool('close')
    matlabpool(nprocs)
end

try
    time=nc_varget(ncfile,'time');
catch
    time=nc_varget(ncfile,'Time');
end

nmonths=length(time);
num_proc_blocks=floor(nmonths)/months_per_proc;
if num_proc_blocks-floor(num_proc_blocks)>0;
    num_proc_blocks=floor(num_proc_blocks)+1;
end
tic;

mean_ssh=zeros(210,360);

parfor block=1:num_proc_blocks   
    disp(sprintf('Block %d/%d',block,num_proc_blocks))
    startidx=(block-1)*months_per_proc;
    start=[ startidx 0 0 0 ];    
    
    if block==num_proc_blocks
        count=[ length(startidx:nmonths)-1 -1 -1 -1 ];
    else
        count=[ months_per_proc -1 -1 -1];
    end
    ssh=squeeze(sum(sum(nc_varget(ncfile,field,start,count),2),1));
    mean_ssh=mean_ssh+ssh;
    
end

mean_ssh=mean_ssh./nmonths;

toc;
matlabpool close