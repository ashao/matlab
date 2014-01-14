function [ ] = run_longterm(numproc,inpath,outpath)

	matlabpool(numproc)
	files=dir([inpath filesep '*.mat']);
	length(files)
	parfor i=1:length(files)
		optcruncher_2011(inpath,outpath,i);
	end
	matlabpool close

end
