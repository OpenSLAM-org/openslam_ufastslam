function particle = predictACFRu(particle,V,G,Q,vehicle,dt,n,lambda,wg,wc)

% INPUTS:
% xv - vehicle pose sample
% Pv - vehicle pose predict covariance
% Note: Pv must be zeroed (or near zero) after each observation. It accumulates the
% vehicle pose uncertainty between measurements.

xv= particle.xv;
Pv= particle.Pv;
dimv= size(xv,1);
dimQ= size(Q,1);

% state augmentation: process noise only
x_aug = [xv ; zeros(dimQ,1)]; 
P_aug = [Pv zeros(dimv,dimQ) ; zeros(dimQ,dimv) Q]; 

% set sigma points
Z = (n+lambda)*(P_aug) + eps*eye(n);
S = chol(Z)';

Kaix = zeros(n,2*n+1);
Kaix(:,1) = x_aug;
for k=1:n
    Kaix(:,k+1) = x_aug + S(:,k);
    Kaix(:,k+1+n) = x_aug - S(:,k);    
end

Kaiy = zeros(size(xv,1),2*n+1);
for k=1:(2*n+1) % transform sigma pts through the process model
    
    Vn = V + Kaix(4,k); % add process noise of linear speed if exist in Kaix
    Gn = G + Kaix(5,k); % add process noise of steering if exist in Kaix
    
    Vc = Vn / (1-tan(Gn)*vehicle.H/vehicle.L); % transform
    
    Kaiy(1,k) = Kaix(1,k) + dt*(Vc*cos(Kaix(3,k)) - Vc/vehicle.L*tan(Gn)*(vehicle.a*sin(Kaix(3,k))+vehicle.b*cos(Kaix(3,k)))); 
    Kaiy(2,k) = Kaix(2,k) + dt*(Vc*sin(Kaix(3,k)) + Vc/vehicle.L*tan(Gn)*(vehicle.a*cos(Kaix(3,k))-vehicle.b*sin(Kaix(3,k))));
    Kaiy(3,k) = Kaix(3,k) + Vc*dt*tan(Gn)/vehicle.L;
end
 
xv_p = 0;
for k=1:(2*n+1)
    xv_p = xv_p + wg(k)*Kaiy(:,k);
end

Pv_p = 0;
for k=1:(2*n+1)
    Pv_p = Pv_p + wc(k)*(Kaiy(:,k) - xv_p)*(Kaiy(:,k) - xv_p)';
end

particle.xv = xv_p;
particle.Pv = Pv_p;
particle.Kaiy = Kaiy; % keep this for a following measurement update