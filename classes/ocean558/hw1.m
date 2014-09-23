%% Model with excess radiation of 2 W/m^2
F.time=(0:499)*86400*365;
F.data=2*ones(size(F.time));
dt=50*86400;
nyears=500; t0=0;

[ Tsurf Tdeep tout ]= hw1_odeint( dt, nyears, t0, F );

clf
subplot(2,2,[1 2])
plot(F.time/86400/365,F.data,'k','LineWidth',2);
title('Radiative Forcing Anomaly');
xlabel('Time'); ylabel('Shortwave Radiation (W/m^2)');

subplot(2,2,3)
plot(tout/86400/365,Tsurf.sol,'k','LineWidth',2);
title('Surface Ocean Temperature Anomaly');
xlabel('Time'); ylabel('Temperature (^\circC)');
ylim([-10 10])
subplot(2,2,4)

plot(tout/86400/365,Tdeep.sol,'k','LineWidth',2);
title('Deep Ocean Temperature Anomaly');
xlabel('Time'); ylabel('Temperature (^\circC)');
ylim([-0.2 1])

%% Model with Gaussian noise to forcing 2 +/- 10 W/m^2
F.time=(0:499)*86400*365;
F.data=2+randn(500,1)*10;
dt=50*86400;
nyears=500; t0=0;

[ Tsurf Tdeep tout ]= hw1_odeint( dt, nyears, t0, F );
clf
subplot(2,2,[1 2])
plot(F.time/86400/365,F.data,'k','LineWidth',2);
title('Radiative Forcing Anomaly');
xlabel('Time'); ylabel('Shortwave Radiation (W/m^2)');
grid on;

subplot(2,2,3)
plot(tout/86400/365,Tsurf.sol,'k','LineWidth',2);
title('Surface Ocean Temperature Anomaly');
xlabel('Time'); ylabel('Temperature (^\circC)');
ylim([-10 10])
grid on;

subplot(2,2,4)
plot(tout/86400/365,Tdeep.sol,'k','LineWidth',2);
title('Deep Ocean Temperature Anomaly');
xlabel('Time'); ylabel('Temperature (^\circC)');
ylim([-0.2 1])
grid on

%% Model with realistic forcing
F20=load('20thcenturyforcingandobs.csv');
data=load('annual.ocean.90S.90N.df_1901-2000mean.dat');
OTemp.time=data(:,1);
OTemp.data=data(:,2);
OTemp.data=OTemp.data-interp1(OTemp.time,OTemp.data,1900);
F.time=F20(:,1)*86400*365;
F.data=F20(:,2);
dt=30*86400;
nyears=length(F.time); t0=F.time(1);

%%
beta=2.0; gamma=1.0;

[ Tsurf Tdeep tout ]= hw1_odeint( dt, nyears, t0, F, beta, gamma,100, 4000 );
fprintf('Beta= %1.2f, Gamma=%1.2f\n',beta,gamma)

clf
subplot(2,2,[1 2])
plot(F.time/86400/365,F.data,'k','LineWidth',2);
title('Radiative Forcing');
xlim([1900 2010])
xlabel('Time'); ylabel('Shortwave Radiation (W/m^2)');
grid on;

subplot(2,2,3); hold on;
plot(tout/86400/365,Tsurf.sol,'k','LineWidth',2);
plot(OTemp.time,OTemp.data,'k--')
title('Surface Ocean Temperature Anomaly');
xlabel('Time'); ylabel('Temperature (^\circC)');
ylim([-1 1]);xlim([1900 2010])
grid on;

subplot(2,2,4)
plot(tout/86400/365,Tdeep.sol,'k','LineWidth',2);
title('Deep Ocean Temperature Anomaly');
xlabel('Time'); ylabel('Temperature (^\circC)');
ylim([-0.1 0.1]);xlim([1900 2010])
grid on

%% Find best fits for beta and gamma