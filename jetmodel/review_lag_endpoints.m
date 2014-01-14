function [  ] = review_lag_endpoints( inpath, start )

files=dir([inpath filesep '*.mat']);
nfiles=length(files);
outpath=[inpath filesep 'confirmed' filesep];
mkdir(outpath);


for ii=start:nfiles
    
    load([inpath filesep files(ii).name]);
    disp(sprintf('File %d',ii))
    if exist('track')
        array=track;
        clear track;
    end
    ndel=0;
    delidx=[];
    
    %% Loop over every entered center
    [ncenters ~]=size(array.optpar_pos);
    disp(sprintf('File %d has %d centers',ii,ncenters))
    for jj=1:ncenters
        %% Positive Events
        subplot(2,1,1,'replace'); hold on;
        plot(array.lat,array.skewness_pos,'k')
        plot(array.optpar_pos(jj,:),[0 1 -1],'-x')
        xlim([-65 -30])
        ylim([-1.5 1.5])
        grid on
        %% Negative Events
        subplot(2,1,2,'replace'); hold on;
        plot(array.lat,array.skewness_neg,'k')
        plot(array.optpar_neg(jj,:),[0 1 -1],'-x')
        xlim([-65 -30])
        ylim([-1.5 1.5])
        grid on
        %% Choose whether to change center
        choice=input('Keep/Delete/Change (1/0/-1)? ');
        if isempty(choice)
            choice=1;
        end
        switch choice
            case 0
                ndel=ndel+1;
                delidx(ndel)=jj;
            case -1
                [centerlat ~]=ginput(2);
                array.optpar_pos(jj,1)=centerlat(1);
                [array.optpar_pos(jj,2) array.optpar_pos(jj,3)]=set_endpoints(array.lat,array.skewness_pos,centerlat(1));
                array.optpar_neg(jj,1)=centerlat(2);
                [array.optpar_neg(jj,2) array.optpar_neg(jj,3)]=set_endpoints(array.lat,array.skewness_neg,centerlat(2));
        end
    end
    
    if ~isempty(delidx)
        array.optpar_pos(delidx,:)=[];
        array.optpar_neg(delidx,:)=[];
    end
    
    if ~isempty(array.optpar_pos)
        save([outpath filesep files(ii).name],'array');
    end
    
    
end

