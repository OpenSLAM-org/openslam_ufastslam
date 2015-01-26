function data = uslam

clear; setconfigfile; h= setup_animations;
particles = initialize_particles(NPARTICLES);
profile off;
profile on -detail builtin

for t=1:n_sampling-1, % for all sampling steps (use faster way)   
    if V(t) ~= 0,
        
        % Predict vehicle state
        for i=1:NPARTICLES
            particles(i)= predictACFRu(particles(i),V(t),G(t),Qe,vehicle,dt,n_r,lambda_r,wg_r,wc_r);
        end
        
        % Measurement update
        for i_obs= featurecount:n_lasersampling,
            
            % Find measurements according to the proper time sequence
            if i_obs < n_lasersampling && lasersampling(i_obs) <= sampling(t+1) && lasersampling(i_obs) > sampling(t),
                
                % Get observation
                z = getObservation(featurecount,OD,OA,perceplimit);
                if ~isempty(z), plines= make_laser_lines(z,particles(1).xv); end % use the first particle for drawing the laser line                               
                
                % Data associations(per-particle)
                for i=1:NPARTICLES                       
                    [particles(i).zf, particles(i).idf, particles(i).zn]= data_associateNN(particles(i), z, Re,GATE_REJECT,GATE_AUGMENT_NN); 
%                     [particles(i).zf, particles(i).idf, particles(i).zn]= data_associate(particles(i), z, Re, GATE_REJECT, GATE_AUGMENT); 
                end
                
                % Known map features
                for i=1:NPARTICLES
                    if ~isempty(particles(i).zf)                      
                        % Sample from optimal proposal distribution
                        particles(i) = sample_proposaluf(particles(i),particles(i).zf,particles(i).idf,Re,n_aug,lambda_aug,wg_aug,wc_aug);                                                
                        % Map update
                        particles(i)= feature_updateu(particles(i),particles(i).zf,particles(i).idf,Re,n_f_a,lambda_f_a,wg_f_a,wc_f_a);                        
                    end
                end 
                
                % Resampling *before* computing proposal permits better particle diversity
                particles= resample_particles(particles, NEFFECTIVE);            
                
                % When new features are observed, augment it to the map
                for i=1:NPARTICLES        
                    if ~isempty(particles(i).zn)
                        if isempty(particles(i).zf) % sample from proposal distribution (if we have not already done so above)
                            particles(i).xv= multivariate_gauss(particles(i).xv, particles(i).Pv, 1);
                            particles(i).Pv= eps*eye(3);
                        end                        
                        particles(i)= add_feature(particles(i), particles(i).zn,Re);
                    end
                end
                
                featurecount = i_obs;            
                break;
            end            
        end
        
        % Plot
        epath = get_epath(particles, epath, NPARTICLES);
        do_plot(h, particles, plines, epath);        
    end
end

profile report

data= particles;

% 
%

function p = initialize_particles(np)
for i=1:np
    p(i).w= 1/np; % initial particle weight
    p(i).xv= zeros(3,1); % initial vehicle pose
    p(i).Pv= eps*eye(3); % initial robot covariance that considers a numerical error
    p(i).Kaiy= []; % temporal keeping for a following measurement update
    p(i).xf= []; % feature mean states
    p(i).Pf= []; % feature covariances    
    p(i).zf= []; % known feature locations
    p(i).idf= []; % known feature index
    p(i).zn= []; % New feature locations   
end

function p= make_laser_lines (rb,xv)
if isempty(rb), p=[]; return, end
len= size(rb,2);
lnes(1,:)= zeros(1,len)+ xv(1);
lnes(2,:)= zeros(1,len)+ xv(2);
lnes(3:4,:)= TransformToGlobal([rb(1,:).*cos(rb(2,:)); rb(1,:).*sin(rb(2,:))], xv);
p= line_plot_conversion (lnes);

function p= make_covariance_ellipses(particle)
N= 10;
inc= 2*pi/N;
phi= 0:inc:2*pi;
circ= 2*[cos(phi); sin(phi)];
p= make_ellipse(particle.xv(1:2), particle.Pv(1:2,1:2) + eye(2)*eps, circ);
lenf= size(particle.xf,2);
if lenf > 0  
    xf= particle.xf;
    Pf= particle.Pf;
    p= [p zeros(2, lenf*(N+2))];
    ctr= N+3;
    for i=1:lenf
        ii= ctr:(ctr+N+1);
        p(:,ii)= make_ellipse(xf(:,i), Pf(:,:,i), circ);
        ctr= ctr+N+2;
    end
end

function p= make_ellipse(x,P,circ)
% make a single 2-D ellipse 
r= sqrtm_2by2(P);
a= r*circ;
p(2,:)= [a(2,:)+x(2) NaN];
p(1,:)= [a(1,:)+x(1) NaN];

function h= setup_animations
figure;
axis([-150 250 -100 250]);
xlabel('[m]'); ylabel('[m]');
hold on, axis equal
h.obs= plot(0,0,'g','erasemode','xor'); % observations
h.xfp= plot(0,0,'r.','erasemode','background'); % estimated features (particle means)
h.xvp= plot(0,0,'r.','erasemode','xor'); % estimated vehicle (particles)
h.cov= plot(0,0,'erasemode','xor'); % covariances of max weight particle
h.epath= plot(0,0,'k','erasemode','xor'); 

function do_plot(h, particles, plines, epath)
xvp = [particles.xv];
xfp = [particles.xf];
w = [particles.w]; 
ii= find(w== max(w)); 
if ~isempty(xvp), set(h.xvp, 'xdata', xvp(1,:), 'ydata', xvp(2,:)), end
if ~isempty(xfp), set(h.xfp, 'xdata', xfp(1,:), 'ydata', xfp(2,:)), end
if ~isempty(plines), set(h.obs, 'xdata', plines(1,:), 'ydata', plines(2,:)), end
pcov= make_covariance_ellipses(particles(ii(1)));
if ~isempty(pcov), set(h.cov, 'xdata', pcov(1,:), 'ydata', pcov(2,:)); end
set(h.epath, 'xdata', epath(1,:), 'ydata', epath(2,:))
drawnow
