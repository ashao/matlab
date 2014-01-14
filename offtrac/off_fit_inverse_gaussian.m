function [ optimal, data ] = off_fit_inverse_gaussian( time, data )
%% Fits an inverse gaussian to the data using constrained simulated annealing

lb=[0 0 0 -Inf];
ub=[Inf Inf Inf Inf];


figure(1)
clf
hold on

iter=0;
itermax=1000;
lambda=[2 200*10 200*10];
tol=inf;
dlambda=inf;
% options=saoptimset('Display','Iter');
% optimal = simulannealbnd(@calc_fit,[1 1 1 5e3],lb,ub,options);
% disp(optimal)

while iter<itermax | tol > 1e-6 | max(dlambda) > 1e-6


alpha=lambda(1);
gamma=lambda(2);
delta=lambda(3);

dbeta=data-inverse_gaussian(lambda,time);
    
Jalpha=sqrt(gamma^3./(4*pi*delta^2*time.^3)).*...
    exp(-gamma*(time-gamma).^2./(4*delta^2*time));
Jgamma=exp(-gamma*(time-gamma).^2./(4*delta^2*time)).*...
    sqrt(gamma).*(time.^2*gamma-4*time*gamma^2+3*gamma^3-6*time*delta^2);
Jdelta=exp(-gamma*(time-gamma).^2./(4*delta^2*time)).*...
    gamma^(3/2).*((time.^2)*gamma+gamma^3-2*time*(gamma^2+delta^2));


A=[Jalpha Jgamma Jdelta];
dlambda=(A'*dbeta)\(A'*A);
lambda=lambda+dlambda;
tol=sum(dbeta);

fprintf('Iter: %d Tol: %f\n',iter, tol)
disp(lambda)
iter=iter+1;

end
