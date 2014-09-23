load /ltraid4/aviso/alongtrack/sla/vxxc_matlab/t001.mat

[ntime npts] = size(track.sla);
ntests = 100;
skew = zeros(ntests,npts);

for i=1:ntests
   
    avgidx = randi(ntime,100,1);
    skew(i,:) = skewness(track.sla(avgidx,:));
    
    
end


clf; hold on
plot(track.lat,std(skew))
plot(track.lat,skewness(track.sla),'--r')
xlim([-65 -30]); ylim([-1 1])