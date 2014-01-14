function [ cfc11 cfc12 sf6 ] = make_CFC_SF6_tables(atmhist)

% atmhist = atmospheric_histories;
gamma = 0.1:0.1:200;
delta = gamma*2;
ngamma=length(gamma);

cfc11.gamma = gamma;
cfc11.delta = delta;

cfc12.gamma = gamma;
cfc12.delta = delta;

sf6.gamma = gamma;
sf6.delta = delta;


t=0:0.01:2000;
measyears = (1980:2010);
measyears = 2010;
cfc11.years = measyears;
cfc12.years = measyears;
sf6.years = measyears;

nyears = length(measyears);

cfc11.ttdconv=zeros(nyears,ngamma);
cfc12.ttdconv=zeros(nyears,ngamma);
sf6.ttdconv=zeros(nyears,ngamma);

for gidx = 1:ngamma
    
    G = inverse_gaussian_waugh([gamma(gidx) delta(gidx)],t);
    
    for yidx = 1:nyears
        
        timesince = measyears(yidx)-atmhist.time;
        maxtime = max(timesince);
        tidx = t<=maxtime & t>=0;
        t_trunc=t(tidx);
        G_trunc=G(tidx);
        
        source.cfc11 = interp1(timesince,atmhist.cfc11.nval,t_trunc);
        source.cfc12 = interp1(timesince,atmhist.cfc12.nval,t_trunc);
        source.sf6 = interp1(timesince,atmhist.sf6.nval,t_trunc);
        
        cfc11.ttdconv(yidx,gidx)=trapz(t_trunc,source.cfc11.*G_trunc);
        cfc12.ttdconv(yidx,gidx)=trapz(t_trunc,source.cfc12.*G_trunc);
        sf6.ttdconv(yidx,gidx)=trapz(t_trunc,source.sf6.*G_trunc);
        
    end
    
        if fix(gamma(gidx))==gamma(gidx)
            fprintf('Year %d Gamma %d CFC-11: %f CFC-12: %f SF6: %f\n', ...
                floor(measyears(yidx)), gamma(gidx), cfc11.ttdconv(yidx,gidx), ...
                cfc12.ttdconv(yidx,gidx),sf6.ttdconv(yidx,gidx));
        end
    
end

end