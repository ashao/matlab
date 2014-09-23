function [ Tsurf Tdeep tout]= hw1_odeint( dt, nyears, t0, F, beta, gamma, h1, h2 )

if nargin<8
    h2=4000;
    if nargin< 7
        h1=50;
        if nargin < 6
            gamma=0.5;
            if nargin <5
                beta=1.2;
            end
        end
    end
end




endt=86400*365*nyears; % End time
nsteps=endt/dt;
disp([h1 h2])
c=h1*1025*3985;
c0=h2*1025*3985;

% Define the dT=f(t)*dt equations
Tsurf.f=inline('1/c*(-beta*Tsurf-gamma*(Tsurf-Tdeep)+F)','Tsurf','Tdeep',...
    'F','c','beta','gamma');
Tdeep.f=inline('gamma/c0*(Tsurf-Tdeep)','Tsurf','Tdeep','gamma','c0');

Tsurf.init=0;
Tdeep.init=0;

% Preallocate solution vectors
Tsurf.sol=zeros(nsteps,1);
Tdeep.sol=zeros(nsteps,1);
tout=zeros(nsteps,1);

% Preallocate coefficient matrix for Runge-Kutta
k=zeros(4,2);

% Set initial conditions
Tsurf.sol(1)=Tsurf.init;
Tdeep.sol(1)=Tdeep.init;

time=t0; tout(1)=t0;
% fprintf('Iter\tTime\tTsurf\tTdeep\n')
for t=1:(nsteps-1)
    
    Fint=interp1(F.time,F.data,time);
    k(1,1)=dt*Tsurf.f(Tsurf.sol(t),Tdeep.sol(t),Fint,c,beta,gamma);
    k(1,2)=dt*Tdeep.f(Tsurf.sol(t),Tdeep.sol(t),gamma,c0);
    
    Fint=interp1(F.time,F.data,time+dt/2);
    k(2,1)=dt*Tsurf.f(Tsurf.sol(t)+0.5*k(1,1),Tdeep.sol(t)+0.5*k(1,2), ...
        Fint,c,beta,gamma);
    k(2,2)=dt*Tdeep.f(Tsurf.sol(t)+0.5*k(1,1),Tdeep.sol(t)+0.5*k(1,2), ...
        gamma,c0);
    
    k(3,1)=dt*Tsurf.f(Tsurf.sol(t)+0.5*k(2,1),Tdeep.sol(t)+0.5*k(2,2), ...
        Fint,c,beta,gamma);
    k(3,2)=dt*Tdeep.f(Tsurf.sol(t)+0.5*k(2,1),Tdeep.sol(t)+0.5*k(2,2), ...
        gamma,c0);
    
    Fint=interp1(F.time,F.data,time+dt);
    k(4,1)=dt*Tsurf.f(Tsurf.sol(t)+k(3,1),Tdeep.sol(t)+k(3,2), ...
        Fint, c, beta, gamma);
    k(4,2)=dt*Tdeep.f(Tsurf.sol(t)+k(3,1),Tdeep.sol(t)+k(3,2), ...
        gamma,c0);
    
    Tsurf.sol(t+1)=Tsurf.sol(t)+(k(1,1)+2*k(2,1)+2*k(3,1)+k(4,1))/6;
    Tdeep.sol(t+1)=Tdeep.sol(t)+(k(1,2)+2*k(2,2)+2*k(3,2)+k(4,2))/6;
    if mod(t,100)==0
        %         fprintf('%d\t%d\t%.2f\t%.2f\n',t+1,floor(time/86400/365),...
        %             Tsurf.sol(t+1),Tdeep.sol(t+1))
    end
    time=time+dt;
    tout(t+1)=time;
end