function [ ] = optcruncher( inpath, outpath, start, finish )



format long

if ~exist(outpath,'dir')
mkdir(outpath)
end
files=dir([inpath '*.mat']);
%     delete([outpath '*.mat']);
%     files
if nargin<3
    start=1;
    finish=length(files);
end

input.number=double(0);
input.optpar=double(0);
input.lat=double(0);
input.skewness=double(0);

for i=start:finish
    
    load([inpath files(i).name]);
    if exist('track')
        array=track;
    end
    %        input.number=double(char(files(i).name(1:4)));
    [m n]=size(array.optpar_pos);
    for j=1:m
        minlat=array.optpar_pos(j,2);
        maxlat=array.optpar_pos(j,3);
        indices=find(array.lat>=minlat & array.lat<=maxlat);
	if isempty(indices) | length(indices)>100
		disp('Skipping')
		break
	end
        input.optpar=double(array.optpar_pos(j,:));
        input.number=files(i).name;
        input.lat=double(array.lat(indices));
        input.skewness=double(array.skewness_pos(indices));
	array.optimal_pos(j,:)=skewfit_ga(input.lat,input.skewness,double(array.optpar_pos(j,1)));
        width=sqrt(2*array.optimal_pos(j,1)+array.optimal_pos(j,3)^2);
        center=array.optimal_pos(j,2);
        SNR=array.optimal_pos(j,4);
        disp(sprintf('Track %d Positive Parameters (Center/Width/SNR): %f / %f / %f',i,center,width,SNR));
        
        minlat=array.optpar_neg(j,2);
        maxlat=array.optpar_neg(j,3);
        indices=find(array.lat>=minlat & array.lat<=maxlat);
	if isempty(indices) | length(indices)>100
		disp('Skipping')
		break
	end
        input.optpar=double(array.optpar_neg(j,:));
        input.number=files(i).name;
        input.lat=double(array.lat(indices));
        input.skewness=double(array.skewness_neg(indices));
        array.optimal_neg(j,:)=skewfit_ga(input.lat,input.skewness,double(array.optpar_neg(j,1)));       
        width=sqrt(2*array.optimal_neg(j,1)+array.optimal_neg(j,3)^2);
        center=array.optimal_neg(j,2);
        SNR=array.optimal_neg(j,4);
        disp(sprintf('Track %d Negative Parameters (Center/Width/SNR): %f / %f / %f',i,center,width,SNR));
    end
    if isfield(array,'optimal_neg')
    save([outpath files(i).name],'array');
	end
end

end


