function [ optimal,fitness,exitflag ] = skewfit_sa( array, numtracks, randn_seed )

    minlat=0;
    maxlat=1;
    LB=[0.001 0 0.001];
    UB=[10 1 10];
    null=[];
    x0=rand(1,3);

    randn_seed=randn(numtracks,1);
    while abs(skewness(randn_seed))>0.001
        randn_seed=randn(numtracks,1);
    end

    simlat=array.lat(array.lat >= minlat & array.lat <=maxlat);
    skewrange=array.skewness(array.lat >= minlat & array.lat<=maxlat);
    options=saoptimset('Display','Final','AnnealingFcn',@annealingboltz,'HybridFcn',@fmincon);
    [optimal residga]=simulannealbnd(@fitskew,x0,LB,UB,options);

    function residual = fitskew(parameters)

        simskew=gensim(parameters);
        residual=sum(sqrt((simskew-skewrange).^2));

    end

    function simskew = gensim(parameters)

        var=parameters(1);
        center_x0=parameters(2);
        L=parameters(3);
        x0=center_x0+sqrt(var)*randn_seed;
        nl=length(randn_seed);
        track_heights=erf((ones(nl,1)*simlat-x0*ones(1,length(simlat)))/(L));
        simskew=skewness(track_heights);

    end


end

function [ randn_seed ] = makerandn( numtracks )
    skewdist=Inf;
    %     while abs(skewdist)>0.01
    randn_seed=randn(numtracks,1);
    skewdist=skewness(randn_seed);
    %     end
end
