function epath = get_epath(particles, epath, NPARTICLES)

% vehicle state estimation result
xvp = [particles.xv];
w = [particles.w]; 
ws= sum(w); w= w/ws; % normalize

% weighted mean vehicle pose
xvmean= 0;
for i=1:NPARTICLES,
    xvmean = xvmean+w(i)*xvp(:,i);               
end

% keep the pose for recovering estimation trajectory
epath=[epath xvmean]; 
