%% BME430 Project 
clc;clear

% Specify parameters
params(1,1) = 1; % kM
params(2,1) = 1; % Erk_Max
params(3,1) = 1; % PI3K_Max
params(4,1) = 360;% kAP, controls divergence of Akt
params(5,1) = .01;% KAP
params(6,1) = 72; % kAD
params(7,1) = 0.01; % KAD
params(8,1) = 360; % kGP
params(9,1) = 0.01; % KGP
params(10,1) = 72; % kGD
params(11,1) = 0.01; % KGD
params(12,1) = 2.3; % kMS, controls Myc jump at start %%%%%%%
params(13,1) = 0.01; % KMS, controls Myc jump at start, change does little
params(14,1) = 0.4; % kMT, at sufficient threshold, increases 58, but does it too early; need help from 62 constants?
params(15,1) = 0.01; % KMT
params(16,1) = 2.08; % dM, correct, controls Myc ss
params(17,1) = 0.35; % dMS, correct
params(18,1) = 2.08; % dMT, correct
params(19,1) = 1;


% params(1,1) = 1; % kM
% params(2,1) = 1; % Erk_Max
% params(3,1) = 1; % PI3K_Max
% params(4,1) = 36;% kAP, controls divergence of Akt
% params(5,1) = .1;% KAP
% params(6,1) = 72e-5; % kAD
% params(7,1) = 0.01e-6; % KAD
% params(8,1) = 36; % kGP
% params(9,1) = 0.01e-6; % KGP
% params(10,1) = 72e-5; % kGD
% params(11,1) = 0.01e-6; % KGD
% params(12,1) = 6.7; % kMS, controls Myc jump at start %%%%%%%
% params(13,1) = 0.1e-2; % KMS, controls Myc jump at start, change does little
% params(14,1) = 0.4e2; % kMT, at sufficient threshold, increases 58, but does it too early; need help from 62 constants?
% params(15,1) = 0.3e-2; % KMT
% params(16,1) = 4.5; % dM, correct, controls Myc ss
% params(17,1) = 1.4; % dMS, correct
% params(18,1) = 5; % dMT, correct
% params(19,1) = 3;

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
% Erk(tspan >= 3 & tspan <= 7) = 1;


% Specify initial values
initvalue = zeros(diffeqn,1); 
initvalue(1,1) = 0.25;
initvalue(2,1) = 0.6;
initvalue(4,1) = 0.6; %initvalue(4,1) = 0.6e-6;
% initvalue(8,1) = 1; %i

% Refer to file with ODEs and call solver
[tsim, results] = ode15s(@BME430Project_core,tspan,initvalue,options,params,PI3K,Erk);
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
