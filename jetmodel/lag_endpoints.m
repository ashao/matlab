function [ ] = lag_endpoints( initpath, lagpath, outpath )

	initfiles=dir([initpath filesep '*.mat']);
	lagfiles=dir([lagpath filesep '*.mat']);
	mkdir(outpath)
	
	for i=1:length(initfiles)
	
		load([initpath filesep initfiles(i).name])
		inittrack=track;clear track;
		load([lagpath filesep lagfiles(i).name])
		track=array;clear array;
		[nguess ~]=size(inittrack.optpar);
		ncenters=0;

		for j=1:nguess
			clf
			subplot(2,1,1);hold on;
			notnan=find(~isnan(track.skewness_pos));
			skew=track.skewness_pos(notnan);
			lat=track.lat(notnan);
			filtskew=filtfilt(ones(5,1),1,skew)/5^2;
			plot(track.lat,track.skewness_pos)
			plot(lat,skew,'LineWidth',2);
			scatter(inittrack.optpar(j,1),0)
			scatter(inittrack.optpar(j,2),-1*ones(size(inittrack.optpar(j,2))),'x')
			scatter(inittrack.optpar(j,3),-1*ones(size(inittrack.optpar(j,3))),'x')
			xlim([-70 -30])
			ylim([-1.5 1.5])
			grid on
			title(sprintf('Negative Events %d (%d/%d)',i,length(initfiles)))
		
			subplot(2,1,2);hold on
			notnan=find(~isnan(track.skewness_neg));
			skew=track.skewness_neg(notnan);
			lat=track.lat(notnan);
			filtskew=filtfilt(ones(5,1),1,skew)/5^2;
			plot(track.lat,track.skewness_neg)
			plot(lat,skew,'LineWidth',2);
			scatter(inittrack.optpar(j,1),0)
			scatter(inittrack.optpar(j,2),-1*ones(size(inittrack.optpar(j,2))),'x')
			scatter(inittrack.optpar(j,3),-1*ones(size(inittrack.optpar(j,3))),'x')
			xlim([-70 -30])
			ylim([-1.5 1.5])
			grid on
			title(sprintf('Negative Events (%d/%d)',i,length(initfiles)))

			[manualcenter ~]=ginput(1);
			if ~isempty(manualcenter)
				centerlat=manualcenter;
			else
				centerlat=inittrack.optpar(j,1);
			end
		
			deletecenter=input('Delete this center? [y/(n)]','s');
			if ~strcmp(deletecenter,'y')

				ncenters=ncenters+1;			
				track.optpar_pos(ncenters,1)=centerlat;
				[track.optpar_pos(ncenters,2) track.optpar_pos(j,3)]=set_endpoints(track.lat,track.skewness_pos,centerlat);
							
				track.optpar_neg(ncenters,1)=centerlat;
				[track.optpar_neg(ncenters,2) track.optpar_neg(j,3)]=set_endpoints(track.lat,track.skewness_neg,centerlat);				
			end
		end	
		save([outpath filesep initfiles(i).name],'track');					

	end

end
