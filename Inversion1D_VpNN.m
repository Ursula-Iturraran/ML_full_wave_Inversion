% Machine learning as seismic velocity-model building method for full-waveform
% inversion: A case study from the Middle Magdalena Valley in Colombia
% Authors: Ursula Iturraran-Viveros, Andres M. Munoz-Garcia and Octavio Castillo-Reyes

% Sample code to perform 1-D inversion using Non-Linear Least Squares and the
% Kennett analytical solution for 1-D layered medium (used as forward modeling tool).
% Since the data shown in the paper is proprietary information we only give
% an example that is not related to this information but it is suitable to
% illustrate how to use this code.

%Load the parameters: Time, Thickness, Ricker-Wavelet, Rho and Vp
data=load('TimeThickRicSeismicRhoVel.dat');

% Set the different variables
time=data(:,1);     % Time vector
dt=time(2)-time(1); % Delta time
thick=data(:,2);    % Thickness vector
ricker=data(:,3);   % Wavelet extracted from the seismic data
seismic=data(:,4);  % Seismic trace
rho=data(:,5);      % rho=density used to compute the synthetic seismogram
Vp=data(:,6);       % Vp=P-wave velocity used to compute the synthetic seismogram
NZ=size(data,1);    % Number of depth-time samples
ntraces=1;          % Number of traces

wlet=10*ricker;     % Scale the wavelet

% Parameters used by the Kennett routine (Written by T. Mukerji)
C(1) = dt;
C(2) = 2;
C(3) = 0;
C(4) = -1;

% Using Kennett routine compute and plot the synthetic seismogram
Wave = zeros; %Vector of the synthetic traces
% The outer for loop is commented since we have only one seismic trace.
% Uncomment if you have more traces and make the adecuate changes to the
% arrays of densities and intial velocites Vp to include densities and
% velocities for each traces

for i=1 : ntraces
    Layer=[Vp rho thick]; % Physical properties of the 1-D layer medium
    Wave = kennetTM(Layer,wlet,C(1),C(2),C(3),C(4)); % 1-D Syhtetic trace
end

% Plot
figure
plot(time, Wave);
xlabel('time');
ylabel('Synthetic trace');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Do the 1-D inversion for the seismic trace  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = zeros;
r = zeros;

options = optimoptions('lsqnonlin','Display','iter','Algorithm','levenberg-marquardt','TolFun',10^-8,'TolX',10^-8,'MaxIter',2000,'MaxFunEvals',10000000);

vmin=2.2663e+03;
vmax=4.3206e+03;
vminVEC=vmin*ones(size(Vp));
vmaxVEC=vmax*ones(size(Vp));

for i=1 : ntraces
    % Observed seismic trace
    WaveObs = seismic;
    % Time domain
    rho(1:30) = 0;
    thick(1:30) = 0;
    wlet(1:30) = 0;
    misfit=@(vel)kennetTM([vel' rho thick],wlet,C(1),C(2),C(3),C(4))-WaveObs';
    % Time domain
    [newV,resnorm,resid,exitflag,output,lambda,J]=lsqnonlin(misfit,Vp',vminVEC,vmaxVEC,options);
end

% Compare solution with the observed seismic trace
LayerNEW=[newV' rho thick]; % Physical properties of the 1-D layer medium
WaveNEW = kennetTM(LayerNEW,wlet,C(1),C(2),C(3),C(4)); % 1-D Syhtetic trace

% Plot
figure
plot(time, WaveNEW);
xlabel('time');
ylabel('Synthetic trace');
hold on;
plot(time,WaveObs,'r');
legend('Synthetic trace','Observed trace');
