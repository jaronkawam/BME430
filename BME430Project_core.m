function dydt = BME430Project_core(t,y,params,PI3K,Erk)

% define state variables
Myc = y(1);
Akt = y(2);
Akt_p = y(3);
Gsk3b = y(4);
Gsk3b_p = y(5);
Myc_Ser62 = y(6);
Myc_Thr58 = y(7);

step = 0.01;
stop = 15;
tspan = 0:step:stop;
PI3K_t = interp1(tspan, PI3K, t, 'linear', 'extrap');
Erk_t = interp1(tspan, Erk, t, 'linear', 'extrap');

% define parameter values
kM = params(1,1);
Erk_Max = params(2,1);
PI3K_Max = params(3,1);
kAP = params(4,1);
KAP = params(5,1);
kAD = params(6,1);
KAD = params(7,1);
kGP = params(8,1);
KGP = params(9,1);
kGD = params(10,1);
KGD = params(11,1);
kMS = params(12,1);
KMS = params(13,1);
kMT = params(14,1);
KMT = params(15,1);
dM = params(16,1);
dMS = params(17,1);
dMT = params(18,1);
GF = params(19,1);

% specify ODEs
dydt(1,1)= kM*GF-(kMS*Erk_t*Myc/(KMS+Myc))-dM*Myc;
dydt(2,1)= -(kAP*PI3K_t*Akt/(KAP+Akt))+(kAD*Akt_p)/(KAD+Akt_p); 
dydt(3,1)= kAP*PI3K_t*Akt/(KAP+Akt)-(kAD*Akt_p)/(KAD+Akt_p);%-kGP*Akt_p*Gsk3b/(KGP+Gsk3b); 
dydt(4,1)= -kGP*Akt_p*Gsk3b/(KGP+Gsk3b)+kGD*Gsk3b_p/(KGD+Gsk3b_p);%-kMT*Gsk3b*Myc_Ser62/(KMT+Myc_Ser62); 
dydt(5,1)= kGP*Akt_p*Gsk3b/(KGP+Gsk3b)-kGD*Gsk3b_p/(KGD+Gsk3b_p);
dydt(6,1)= kMS*Erk_t*Myc/(KMS+Myc)-(1-PI3K_t)*(kMT*Gsk3b*Myc_Ser62/(KMT+Myc_Ser62))-dMS*Myc_Ser62; 
dydt(7,1)= (1-PI3K_t)*kMT*Gsk3b*Myc_Ser62/(KMT+Myc_Ser62)-dMT*Myc_Thr58; 

return
