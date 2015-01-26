function [zf,idf,zn]= data_associateNN(particle,z,R, gate1, gate2)
% Simple gated nearest-neighbour data-association. 
% Approximately 1.6 times faster than data_associate.m under a condition of
% the same number of feature.
% Chanki Kim, 2010.

zf= []; zn= []; idf= []; 
      
Nf= size(particle.xf,2); % number of Column. = number of features already in map        
zp= zeros(2,1);

xv= particle.xv;
% linear search for nearest-neighbour, no clever tricks (like a quick
% bounding-box threshold to remove distant features; or, better yet,
% a balanced k-d tree lookup). TODO: implement clever tricks.
for i=1:size(z,2)
    jbest= 0;
    outer= inf;
    
    if Nf~=0,
        % search for spatially (deometrically) nearest neighbour (CK)
        dmin = inf;
        jbest_s= 0;     
        
        for j=1:Nf,
            dx= particle.xf(1,j)-xv(1);
            dy= particle.xf(2,j)-xv(2);
            d2= dx^2 + dy^2;
            d= sqrt(d2);
            zp(1,1)= d; 
            zp(2,1)= pi_to_pi(atan2(dy,dx) - xv(3)); 
            v= z(:,i)-zp; v(2) = pi_to_pi(v(2));
            d2= v'*v;
            if d2 < dmin,
                dmin= d2;
                jbest_s= j;
            end
        end    
        
        % Mahalanobis test for the candidate neighbour (CK)
        nis= compute_association_nis(particle,z(:,i),R, jbest_s); % Nearest Neibour 
        if nis < gate1 % if within gate, store nearest-neighbour
            jbest= jbest_s;
        elseif nis < outer % else store best nis value
            outer= nis;    
        end
    end
    
    % add nearest-neighbour to association list
    if jbest ~= 0
        zf=  [zf  z(:,i)];
        idf= [idf jbest];
    elseif outer > gate2 % z too far to associate, but far enough to be a new feature
        zn= [zn z(:,i)];
    end
end

function nis= compute_association_nis(particle,z,R,idf)
% return normalised innovation squared (ie, Mahalanobis distance) and normalised distance
[z_predict,~,~,Sf,Sv]= compute_jacobians(particle, idf, R);
v= z-z_predict;   % innovation
v(2)= pi_to_pi(v(2));

nis= v'/Sf*v;




