ntest = 100;
est_params = zeros(ntest,3);
params = zeros(ntest,4);

simy = linspace(-500,500,100);
in_params = [100 0 100 0.8];

for test = 1:ntest
    
    
    
    fprintf('Test %d/%d\n',test,ntest)
    track_heights = simulate_SSH(in_params,simy);
    simskew = skewness(track_heights);
    fitidx = abs(simy)<120;    
    x0 = [80 20 60 0.9];
    params(test,:) = jet_model_estimate_curvefit(simy(fitidx),simskew(fitidx),[20 -50 20 0.2],[200 50 200 1],x0);
    
    
end
hist(params(:,2))
nanmean(params)
%%
clf; hold on;
simy = linspace(-500,500,50);


in_params = [100 0 100 0.7];

track_heights = simulate_SSH(in_params,simy);
simskew = skewness(track_heights);
% jet_model_estimate_curvefit(simy, simskew, [0.01 -2 ], [10 2], [1 1])
plot(simy,simskew)
fwhm(simy,cumsum(simskew))
sqrt(in_params(1)^2 + in_params(3)^2)

%% Check dependence on double jet


clf; hold on;
simy = linspace(-500,500,50);


in_params1 = [50 150 50 0.8];
in_params2 = [50 -50 50 0.8];
track_heights1 = simulate_SSH(in_params1,simy);
track_heights2 = simulate_SSH(in_params2,simy);
track_heights = track_heights1 + track_heights2;

simskew = skewness(track_heights);
% jet_model_estimate_curvefit(simy, simskew, [0.01 -2 ], [10 2], [1 1])
plot(simy,simskew,'k','LineWidth',2)
plot(simy,skewness(track_heights1),'r')
plot(simy,skewness(track_heights2),'b')
plot(simy,skewness(track_heights2)+skewness(track_heights1),'g')
fwhm(simy,cumsum(simskew))

%% Check dependence on each of hte four parameters

ntest = 100;
params = zeros(ntest,3);

simy = linspace(-500,500,100);
in_params = [100 0 100 0.6];
options = optimset('Display','off');
SNR = (0:2:100);
var = 20:5:300;
L = 20:5:300;
maxskew = length(var);

for test = length(var):-1:1
    fprintf('Test %d/%d\n',test,length(var))
    maxskew(test) = 0;
    in_params(1) = var(test);
    for i=1:20
       
       track_heights = simulate_SSH(in_params,simy);
       simskew = skewness(track_heights);
       maxskew(test) = maxskew(test) + max(simskew);
        
    end    
    
    
end
maxskew_var = maxskew/10;

in_params = [100 0 100 0.6];
for test = length(L):-1:1
    fprintf('Test %d/%d\n',test,length(var))
    maxskew(test) = 0;
    in_params(3) = L(test);
    for i=1:20
       
       track_heights = simulate_SSH(in_params,simy);
       simskew = skewness(track_heights);
       maxskew(test) = maxskew(test) + max(simskew);
        
    end    
    
    
end

maxskew_L = maxskew/10;
%%

in_params = [100 0 100 0.6];
clear maxskew
for test = length(SNR):-1:1
    fprintf('Test %d/%d\n',test,length(var))
    maxskew(test) = 0;
    in_params(4) = SNR(test)/100;
    for i=1:20
       
       track_heights = simulate_SSH(in_params,simy);
       simskew = skewness(track_heights);
       maxskew(test) = maxskew(test) + max(simskew);
        
    end    
    
    
end

maxskew_SNR = maxskew/10;
%%
subplot(1,2,1); hold on;
plot(L,maxskew_L,'k--','LineWidth',2)
plot(var,maxskew_var,'k','LineWidth',2)
xlabel('\sigma, L [km]')
ylabel('Skewness Magnitude')
axis square
legend('L','\sigma')
grid on
title('(a)')

subplot(1,2,2);
plot(SNR/100,maxskew_SNR,'k','LineWidth',2)
xlim([0 0.9])
xlabel('SNR')
ylabel('Skewness Magnitude')
axis square
grid on
title('(b)')
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/safaltimetry/skewness_parameters.m')