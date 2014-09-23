function [ optimal fval simskew ] = ...
    jet_model_estimate(datay,dataskew,SNR,LB,UB,x0);
% JET_MODEL_ESTIMATE: Estimates center, width, and variability of a
% meandering Gaussian jet based on the skewness fields
%   Algorithm description:
%       A satellite pass of a meandering jet will observe a change in SSH
%       (Gille 1994) of the form
%           SSH(x,t) = A*erf( (x-x0)/L^2)
%       where A is the total height difference across the jet and x0 is a
%       normal random variabile which describes a changing center of the
%       jet, and L is the width of the jet (in statistical terms, the
%       standard deviation of the gaussian-shaped jet
%
%       Per Thompson and Demirov (2005), skewness of sea surface height  is a
%       robust statistic that can be used to identify frontal features.
%       For a eastward jet, skewness should be positive south of the jet
%       and negative north of it.
%
%       Here we use these two results and the skewness calculated from the
%       TOPEX/POSEIDON, Jason-1, and Jason-2 satellite altimetry missions
%       to estimate the position of fronts in the Southern Ocean.

%% Declare some parameters
numtracks = 1000;
% LB=[0.001 -5 0.001];
% UB=[10 5  10];
% x0 = [1.5 0.2 2];

% datap = polyfit(datay,dataskew,1);
% dataline = polyval(datap,datay);


%% Make a random seed with small skewness
randn_seed=randn(numtracks,1);
while abs(skewness(randn_seed))>0.001
    randn_seed=randn(numtracks,1);
end


options=psoptimset('Display','None','TolMesh',1e-12,'TolX',1e-12,'TolFun',1e-9,'CompletePoll','on','CompleteSearch','on');
%  options = saoptimset('Display','Iter','TolFun',1e-9);
% [optimal fval]=simulannealbnd(@fitskew,x0,LB,UB,options);
[optimal fval]=patternsearch(@fitskew,x0,[],[],[],[],LB,UB,[],options);
LB = optimal-optimal*0.1;
UB = optimal+optimal*0.1;
% [optimal fval]=fmincon(@fitskew,optimal,[],[],[],[],LB,UB,[],options);
% simskew = gensim(optimal);

    function residual = fitskew(parameters)

        simskew=gensim(parameters);        
%         simp = polyfit(datay,simskew,1);
%         simline = polyval(simp,datay);
%         residual = sum( (simline-dataline).^2);
        residual = nansum( (dataskew-simskew).^2 );
        
    end

    function simskew = gensim(parameters)
        
        var=parameters(1);
        center_x0=parameters(2);
        L=parameters(3);

        randn_seed=randn(numtracks,1);
        center_x0=center_x0+var*randn_seed;
        nl=length(randn_seed);
%         whos
        track_heights=erf((ones(nl,1)*datay-center_x0*ones(1,length(datay)))/(L));
        track_heights = track_heights + (1-SNR)*randn(size(track_heights));
        simskew=skewness(track_heights);

    end

end
