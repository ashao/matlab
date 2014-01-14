function simskew = gensim(parameters,simlat,randn_heights,randn_x0)
var=parameters(1);
center_x0=parameters(2);
L=parameters(3);
Amp=parameters(4);
x0=center_x0+sqrt(var)*randn_x0;
nl=length(randn_x0);
track_heights=erf((ones(nl,1)*simlat-x0*ones(1,length(simlat)))/(L));
track_heights=Amp*track_heights+randn_heights;
simskew=skewness(track_heights);
end