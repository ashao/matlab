function [ optimal,fitness ] = skewfit_sa( simlat, datskew, center )

warning off all
fminopt=optimset('Display','off','MaxFunEvals',5000);
gaopt=gaoptimset('Display','iter','HybridFcn',{@fmincon,fminopt});
saopt=saoptimset('Display','off','HybridFcn',{@fmincon,fminopt});
psopt=psoptimset('Display','iter','MaxIter',10000,'MaxFunEvals',20000);
LB=[0.001 center-5 0.001 10];
UB=[5 center+5 5 100];

x0=[rand(1) center rand(1) 50];

numtracks=100;
randn_heights=make_zeroskewdist(numtracks,length(simlat));
randn_x0=make_zeroskewdist(1,numtracks)';

% datskew=filtfilt(ones(3,1),1,datskew)/9;
% clf
% hold on
% plot(simlat,datskew)
% [optimal fitness]=ga(@fitskew,4,[],[],[],[],LB,UB,[],gaopt);

for i=1:5
[optimal(i,:) fitness]=simulannealbnd(@fitskew,x0,LB,UB,saopt);
end
% [optimal fitness]=patternsearch(@fitskew,x0,[],[],[],[],LB,UB,[],psopt);
optimal=mean(optimal);

%if abs(optimal(2)-center) > 2
%disp('Checking to validate center')
%[optimal fitness]=simulannealbnd(@fitskew,x0,LB,UB,saopt);
%end

    function residual = fitskew(parameters)
        simskew=gensim(parameters,simlat,randn_heights,randn_x0);
        if isempty(simskew)
            residual=NaN;            
        else
            residual=sum( (simskew-datskew).^2);
        end

    end    
end
