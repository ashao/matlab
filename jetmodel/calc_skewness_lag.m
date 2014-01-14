function [  ] = calc_skewness_lag( allarray, indexval, indexdate, outpath, cutoff )
%% Calculates skewness from climatology index


% +/- 3 month lag in 5 day increments. Positive means that SLA measurement
% is before the index.

lags=-60:5:60; %Setup 5 day leads/lags for up to three months.
mkdir(outpath);
for i=1:length(lags)
    lagdates=indexdate+lags(i); %Add lags
%     sla_index_val=interp1(single(indexdate),single(indexval),allarray.time);
    
    fullpath=[outpath sprintf('lag_%02d',lags(i))];
    mkdir(fullpath);
    
    
    for track=1:254
    disp(sprintf('Calculating skewness for track %d with %d day lag',track,lags(i)))
        array.lat=allarray.latarray(track,:);
        array.lon=allarray.lonarray(track,:);
        array.time=squeeze(allarray.time(:,track,:));
        array.index_lag=lags(i);
        
        track_sla_all=squeeze(allarray.sla(:,track,:));        
        sla_index_val=interp1(lagdates,indexval,array.time);        
	nrecs=sum(~isnan(track_sla_all));
	nrecs(nrecs==0)=NaN;
	nrecs=nanmean(nrecs);	

	ntot=length(find(~isnan(track_sla_all)));        
        
        track_sla_pos=track_sla_all;
        track_sla_pos(sla_index_val<cutoff)=NaN;
	npos=sum(~isnan(track_sla_pos));
	npos(npos==0)=NaN;
	npos=floor(nanmean(npos));	

        track_sla_neg=track_sla_all;
        track_sla_neg(sla_index_val>-cutoff)=NaN;
	nneg=sum(~isnan(track_sla_neg));
	nneg(nneg==0)=NaN;
	nneg=floor(nanmean(nneg));	
        
	disp(sprintf('%d Positive Records / %d Negative Records\n',npos,nneg))
        array.skewness_pos=skewness(track_sla_pos,0,1);
        array.skewness_neg=skewness(track_sla_neg,0,1);
	array=rmfield(array,'time');

        save([fullpath filesep sprintf('t%03d.mat',track)],'array','-v7.3');
    end
        
       
    
    
    

end

