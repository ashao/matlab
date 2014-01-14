function [ G ] = inverse_gaussian( coeffs, time )

mu=coeffs(1);
lambda=coeffs(2);

a1=sqrt(lambda./(2*pi*time.^3));
a2= - lambda*(time-mu).^2./(2*mu^2.*time);

G=a1.*exp(a2);

G(time==0)=0;
