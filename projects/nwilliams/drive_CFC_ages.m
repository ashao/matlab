botdata = load('CFCdataforAndrew.mat');
linenames = fields(botdata);
nlines = length(linenames);
atmhist = atmospheric_histories;
tracprop = tracer_properties;
convK = 273.15;
maxCFC = max(atmhist.cfc11.nval);

trac.atmconc = atmhist.cfc11.nval';
trac.time = atmhist.time;

for iline = [1 5]
    cline = linenames{iline};
    nbot = length(botdata.(cline).SALNTY);
    fprintf('Processing %d bottles from line %s\n', nbot, cline);
    sol = trac_calcsol(botdata.(cline).POTTMP+convK,botdata.(cline).SALNTY, ...
        tracprop.cfc11.sol_coeffs.grav);
    pcfc=botdata.(cline).CFC11./sol;
    pcfc(pcfc>maxCFC)=NaN;
    botdata.(cline).PCFC11=pcfc;
    gamma = zeros(nbot,1);
    ttderr = zeros(nbot,1);
    isvalid = ~isnan(pcfc) & ~isinf(pcfc);
    date = num2str(botdata.(cline).DATE);
    
    parfor ibot = 1:nbot
        
        if isvalid(ibot)
            year = str2double(date(ibot,1:4));
            mon = str2double(date(ibot,5:6));
            day = str2double(date(ibot,7:8));
            meastime = year + (datenum(1,mon,day) - datenum(0,12,31))/365;            
            [mindiff minidx] = min(abs(trac.atmconc-pcfc(ibot)));
            x0 = meastime - trac.time(minidx);
            
            [gamma(ibot) mull ttderr(ibot)] = match1dTTD(...
                pcfc(ibot),0.8,trac,meastime,'x0',x0);
            
            fprintf('Bottle %d/%d pCFC: %f Gamma: %f Error: %e\n', ...
                ibot, nbot, pcfc(ibot),gamma(ibot),ttderr(ibot))
        else
            gamma(ibot)=NaN;
        end

    end
    botdata.(cline).TTD11=gamma;
    botdata.(cline).TTDERR=ttderr;        
    
end
%%
save('sp4.sat80.ttd11.mat','botdata');
%%
save('plines.sat80.ttd11.mat','botdata')

%%
subplot(1,2,1)

scatter(botdata.S4P_1992.TTD11,botdata.S4P_1992.CTDPRS,[],'b','filled')
axis square
grid on
xlabel('TTD Mean Age (CFC-11)')
ylabel('Pressure')
xlim([0 600])
title(linenames(1))
set(gca,'ydir','reverse')

subplot(1,2,2)
scatter(botdata.S4P_2011.TTD11,botdata.S4P_2011.CTDPRS,[],'b','filled')
axis square
grid on
xlim([0 600])
xlabel('TTD Mean Age (CFC-11)')
ylabel('Pressure')
title(linenames(5))
set(gca,'ydir','reverse')