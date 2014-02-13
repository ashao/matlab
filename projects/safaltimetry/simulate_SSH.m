    function [simy track_heights simskew argument] = simulate_SSH(parameters)
        
        var=parameters(1);
        center_y0=parameters(2);
        L=parameters(3);
        
        simy = linspace(center_y0-5*sqrt(L^2+var),center_y0+5*sqrt(L^2+var),100);
        randn_seed=randn(1000,1);
        y0=center_y0+var*randn_seed;
        nl=length(randn_seed);
        size(ones(nl,1)*simy)
        size(center_y0*ones(1,length(simy)))
        argument = (ones(nl,1)*simy-y0*ones(1,length(simy)))/(L);
        track_heights=erf((ones(nl,1)*simy-y0*ones(1,length(simy)))/(L));

%         track_heights=erf(argument.^2);
        track_heights=track_heights+0.2*randn(size(track_heights));
        simskew=skewness(track_heights);

    end