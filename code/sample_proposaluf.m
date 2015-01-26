function particle= sample_proposaluf(particle,zf,idf,R,n,lambda,wg,wc)
% Compute proposal distribution and then sample from it

xv_= particle.xv; % predictive vehicle state
Pv_= particle.Pv; % predictive vehicle state error covariance

lenidf= length(idf); % number of currently observed features
dimv= size(xv_,1); % vehicle state dimension
dimf= size(zf,1); % feature state dimension

% for batch update
z_hat= zeros(dimf*lenidf,1); % predictive observation
z= zeros(dimf*lenidf,1); % sensory observation
A= zeros(dimf*lenidf,2*n+1); % stack of innovation covariance for vehicle uncertainty
wc_s= sqrt(wc);

for i=1:lenidf,
    j=idf(i); % index of this observed feature
    xfi=particle.xf(:,j);  % get j-th feature mean
    Pfi=particle.Pf(:,:,j); % get j-th feature cov.
    z(2*i-1:2*i,1)=zf(:,i); % stack of sensory observations

    % state augmentation
    x_aug = [xv_ ; xfi];
    P_aug = [Pv_ zeros(dimv,dimf) ; zeros(dimf,dimv) Pfi]; 
   
    % set sigma points
    Ps = (n+lambda)*(P_aug) + eps*eye(n); 
    Ss = chol(Ps)';
    Ksi = zeros(n,2*n+1);
    Ksi(:,1) = x_aug;
    for k=1:n
        Ksi(:,k+1) = x_aug + Ss(:,k);
        Ksi(:,k+1+n) = x_aug - Ss(:,k);    
    end
    % passing through observation model
    Ai = zeros(dimf,2*n+1);
    bs=zeros(1,2*n+1); % bearing sign    
    z_hati = 0; % predicted observation ('dimf' by 1)    
    for k=1:(2*n+1) % pass the sigma pts through the observation model
        d= Ksi(dimv+1:dimv+dimf,k) - Ksi(1:dimv-1,k);
        r= sqrt(d(1)^2 + d(2)^2); % range 
        bearing= atan2(d(2),d(1));
        bs(k)=sign(bearing);    
        if k>1, % unify the sign
            if bs(k) ~= bs(k-1);            
                if bs(k)<0 && -pi < bearing && bearing < -pi/2, 
                    bearing=bearing+2*pi;
                    bs(k)=sign(bearing);    
                elseif bs(k)>0 && pi/2 < bearing && bearing < pi,
                        bearing=bearing-2*pi;
                        bs(k)=sign(bearing);    
                end
            end
        end
        Ai(:,k)= [r;      bearing - Ksi(dimv,k)]; % bearing  **do not use pi_to_pi here** 
        z_hati = z_hati + wg(k)*Ai(:,k);  % predictive observation         
    end    
    z_hati_rep= repmat(z_hati,1,2*n+1);
    A(2*i-1:2*i,:)= (Ai - z_hati_rep); 
    for k=1:(2*n+1)
        A(2*i-1:2*i,k) = A(2*i-1:2*i,k) * wc_s(k);
    end    
    z_hati(2)= pi_to_pi(z_hati(2)); % now use pi_to_pi, but this is not imperative.
    z_hat(2*i-1:2*i,1)= z_hati;  
end
% augmented noise matrix
R_aug = zeros(dimf*lenidf,dimf*lenidf);
for i=1:lenidf,
    R_aug(2*i-1:2*i,2*i-1:2*i)= R;
end

% innovation covariance (ISSUE)
S=A*A'; % vehicle uncertainty + map + measurement noise
S=(S+S')*0.5 + R_aug;  % make symmetric for better numerical stability

% cross covariance: considering vehicle uncertainty
X = zeros(dimv,2*n+1); % stack
for k=1:(2*n+1)
    X(:,k)= wc_s(k) * (Ksi(1:3,k) - xv_);     
end

U= X*A'; % cross covariance ('dimv' by 'dimf*lenidf')

% Kalman gain
K = U / S;    

% innovation ('dimf*lenidf' by 1) 
v = z-z_hat;
for i=1:lenidf, 
    v(2*i)=pi_to_pi(v(2*i));
end

% standard Kalman update
xv = xv_ + K*v;
Pv = Pv_ - K*U'; 

% compute weight (parallel process) : ERB for SLAM problem
Lt= S; % square matrix of 'dimf*lenidf'
den= sqrt(2*pi*det(Lt));
num= exp(-0.5 * v' / Lt * v);
w = num/den;
particle.w = particle.w * w;

% sample from proposal distribution
xvs = multivariate_gauss(xv,Pv,1); 
particle.xv = xvs;
particle.Pv= eye(3)*eps; % initialize covariance 

%
%

function x= pi_to_pi(x)
for i=1:size(x,2),
    if x(i) > pi,
        while x(i) > pi
            x(i)= x(i) - 2*pi;
        end
    elseif x(i) < -pi
        while x(i) < -pi
            x(i)= x(i) + 2*pi;
        end
    end
end
