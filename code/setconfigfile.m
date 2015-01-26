%%% Set environment
format compact;
path(path, '../')
plines= zeros(2,1); 
featurecount= 1; % initialize laser feature count
epath= [];  % estimation path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Experiment on Victoria park
load gps_x.txt;
load gps_y.txt;
%%% Entire dataset
load ObservationDistance.txt; OD= ObservationDistance;
load ObservationAngle.txt; OA= ObservationAngle;
load lasersampling.txt;
load time.txt; sampling= time;
load speed.txt; V= speed;
load steering.txt; G= steering;
n_lasersampling = length(lasersampling);
n_sampling = round(length(sampling));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% H/W parameter(truck). Keep the values.
dt= 0.025; % [s], time interval between control signals
vehicle.L= 2.83; % [m]
vehicle.H= 0.76; % [m]
vehicle.b= 0.5; % [m]
vehicle.a= 3.78; % [m]
veh= [0 -vehicle.H -vehicle.H; 0 -1 1];

% control noises. You can change these values for your application.
sigmaV= 2; % [m/s]
sigmaG= (6*pi/180); % [rad]
Qe = [sigmaV^2 0; 0 sigmaG^2];

% observation(measurement) noises 
sigmaR= 1; % [m]
sigmaB= (3*pi/180); % [rad]
Re = [sigmaR^2 0; 0 sigmaB^2];
perceplimit=30; % [m] 

% resampling criteria
NPARTICLES= 5; % number of particles(samples, hypotheses) 
NEFFECTIVE= 0.5*NPARTICLES; % minimum number of effective particles

% data association - innovation gates (Mahalanobis distance)
GATE_REJECT= 5.991; % maximum distance for association
GATE_AUGMENT_NN= 2000; % minimum distance for creation of new feature
GATE_AUGMENT= 100; % minimum distance for creation of new feature (100)

% parameters related to SUT
dimv= 3; dimQ= 2; dimR= 2; dimf= 2;

% vehicle update
n_aug=dimv+dimf; 
alpha_aug=0.9; beta_aug=2; kappa_aug=0;
lambda_aug = alpha_aug^2 * (n_aug + kappa_aug) - n_aug; 
lambda_aug=lambda_aug+dimR;
wg_aug = zeros(1,2*n_aug+1); wc_aug = zeros(1,2*n_aug+1);
wg_aug(1) = lambda_aug/(n_aug+lambda_aug);
wc_aug(1) = lambda_aug/(n_aug+lambda_aug)+(1-alpha_aug^2+beta_aug);
for i=2:(2*n_aug+1)
    wg_aug(i) = 1/(2*(n_aug+lambda_aug));    
    wc_aug(i) = wg_aug(i);    
end

% vehicle prediction 
n_r=dimv+dimQ;
alpha_r=0.9;
beta_r=2;
kappa_r=0;
lambda_r = alpha_r^2 * (n_r + kappa_r) - n_r; 
lambda_r= lambda_r+dimR; % should consider dimension of related terms for obtaining equivalent effect with full augmentation
wg_r = zeros(1,2*n_r+1); wc_r = zeros(1,2*n_r+1);
wg_r(1) = lambda_r/(n_r + lambda_r);
wc_r(1) = lambda_r / (n_r+lambda_r) + (1 - alpha_r^2 + beta_r);
for i=2:(2*n_r+1)
    wg_r(i) = 1/(2*(n_r+lambda_r));    
    wc_r(i) = wg_r(i);    
end

%feature updates (augmented state)
n_f_a= dimf + dimR;  
alpha_f_a=0.9; 
beta_f_a=2; 
kappa_f_a=0;
lambda_f_a = alpha_f_a^2 * (n_f_a + kappa_f_a) - n_f_a;
wg_f_a = zeros(1,2*n_f_a+1); wc_f_a = zeros(1,2*n_f_a+1);
wg_f_a(1) = lambda_f_a / (n_f_a + lambda_f_a);
wc_f_a(1) = lambda_f_a / (n_f_a + lambda_f_a) + (1 - alpha_f_a^2 + beta_f_a);
for i=2:(2*n_f_a+1)
    wg_f_a(i) = 1/(2*(n_f_a + lambda_f_a));    
    wc_f_a(i) = wg_f_a(i);    
end
%feature initialization
n_f= dimf;
alpha_f=0.9; 
beta_f=2; 
kappa_f=0;
lambda_f = alpha_f^2 * (n_f + kappa_f) - n_f;
wg_f = zeros(1,2*n_f+1); wc_f = zeros(1,2*n_f+1);
wg_f(1) = lambda_f / (n_f + lambda_f);
wc_f(1) = lambda_f / (n_f + lambda_f) + (1 - alpha_f^2 + beta_f);
for i=2:(2*n_f+1)
    wg_f(i) = 1/(2*(n_f+lambda_f));    
    wc_f(i) = wg_f(i);    
end
