%% BME430 Project, Fall 2024 
% Partners: C Jiang, M Montiel
% Reconstructed based on Lee et al. 2008
clc;clear;close all

% Specify parameters
params(1,1) = 1; % kM
params(2,1) = 1; % Erk_Max
params(3,1) = 1; % PI3K_Max
params(4,1) = 360;% kAP
params(5,1) = .01;% KAP
params(6,1) = 72; % kAD
params(7,1) = 0.01; % KAD
params(8,1) = 360; % kGP
params(9,1) = 0.01; % KGP
params(10,1) = 72; % kGD
params(11,1) = 0.01; % KGD
params(12,1) = 2.3; % kMS
params(13,1) = 0.01; % KMS
params(14,1) = 0.4; % kMT
params(15,1) = 0.01; % KMT
params(16,1) = 2.08; % dM
params(17,1) = 0.35; % dMS
params(18,1) = 2.08; % dMT
params(19,1) = 1; % GF

% specify solver details
step = 0.01;
stop = 15;
tspan = 0:step:stop;
diffeqn = 7;
options = odeset('RelTol',1e-6);

% Initialize PI3K and Erk
Erk_R = 0.1*params(2,1);
PI3K_R = 0.1*params(3,1);
PI3K = PI3K_R*ones(length(tspan),1);
Erk = Erk_R*ones(length(tspan),1);

% Define PI3K and Erk levels
PI3K(tspan >= 0 & tspan <= 1) = 1;   % PI3K is 1 from 0 to 1
PI3K(tspan >= 4 & tspan <= 7) = 1; % PI3K is 1 again from 3 to 6
Erk(tspan >= 0 & tspan <= 1) = 1; % Erk is 1 from 0 to 1

% Specify initial values
initvalue = zeros(diffeqn,1); 
initvalue(1,1) = 0.25;
initvalue(2,1) = 0.6;
initvalue(4,1) = 0.6;

% Refer to file with ODEs and call solver
[tsim, results] = ode15s(@lee_core,tspan,initvalue,options,params,PI3K,Erk);
Myc_total = results(:,1)+results(:,6)+results(:,7);

% plot results
figure(1);
subplot(3,1,1)
plot(tsim,Erk,'Linewidth',2)
title('Erk')

subplot(3,1,2)
plot(tsim,PI3K,'Linewidth',2)
title('PI3K')

subplot(3,1,3)
plot(tsim,results(:,1),tsim,results(:,6),tsim,results(:,7),tsim,Myc_total,'Linewidth',2)
legend('Myc','Myc Ser62','Myc Thr58','Myc Total')
title('Myc Accumulation')
ylabel('concentration')
xlabel('time')

figure(2)
plot(tsim,results,'Linewidth',2)
legend('Myc','Akt','Akt_p','Gsk3b','Gsk3b_p','Myc Ser62','Myc Thr58')
