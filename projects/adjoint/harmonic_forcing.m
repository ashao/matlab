dt = 0.01;
t = 0:dt:100;
params1 = [ 10 5 ];
params2 = [ 30 10 ];

variability = zeros(size(t));
for i = 1:length(t)

    
    G1 = inverse_gaussian_waugh(params1,t);
    G2 = inverse_gaussian_waugh(params2,t);
    G = 0.5*(G1+G2);
    forcing = cos(2*pi/10 * (t + dt*i));
    variability(i) = trapz(t,G.*forcing);
    
end


% plot(t,forcing)

plot(t,variability);
corrcoef(variability,cos(2*pi/10*t))
%%
varinv = size(variability);
for i = 1:length(t)
   
    G1 = inverse_gaussian_waugh(params1,t);
    G2 = inverse_gaussian_waugh(params2,t);
    G = 0.5*(G1+G2);
    forcing = cos(2*pi/10 * (t - dt*i));
    varinv(i) = trapz(t,fliplr(G).*forcing);
    
end
plot(varinv)