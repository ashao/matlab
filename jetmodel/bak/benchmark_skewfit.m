inpath=sprintf('/scratch/data/jetmodel/toprocess/aao/lag_00/');
outpath='/scratch/data/jetmodel/output/temp/';

maxNumCompThreads(1);
tic
for i=1:5
    optcruncher_2011_diff(inpath,outpath,1,1)
end
disp(sprintf('Average time with one thread: %f',toc/5))

maxNumCompThreads(4);
tic
for i=1:5
    optcruncher_2011_diff(inpath,outpath,1,1)
end
disp(sprintf('Average time with 4 threads: %f',toc/5))
