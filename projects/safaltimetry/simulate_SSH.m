    function [track_heights simskew ] = simulate_SSH(parameters,simy)
        
        var=parameters(1);
        center_y0=parameters(2);
        L=parameters(3);
        SNR = parameters(4);
        numtracks = 1000;
%         simy = linspace(center_y0-5*sqrt(L^2+var),center_y0+5*sqrt(L^2+var),100);
        randn_seed=randn(numtracks,1);
        while abs(skewness(randn_seed))>0.001
            randn_seed=randn(numtracks,1);
        end
        y0=center_y0+var*randn_seed;
        nl=length(randn_seed);       
        
        track_heights=erf((ones(nl,1)*simy-y0*ones(1,length(simy)))/(sqrt(2)*L));

%         track_heights=erf(argument.^2);
        track_heights=SNR*track_heights+(1-SNR)*randn(size(track_heights));
        simskew=skewness(track_heights);

    end