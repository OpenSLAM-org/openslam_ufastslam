function [zf,idf,zn]= data_associate(particle,z,R, gate1, gate2)
% 
% Simple gated nearest-neighbour data-association. 

zf= []; zn= []; idf= []; 
      
Nf= size(particle.xf,2); % number of Column. = number of features already in map        

% linear search for nearest-neighbour, no clever tricks (like a quick
% bounding-box threshold to remove distant features; or, better yet,
% a balanced k-d tree lookup). TODO: implement clever tricks.
for i=1:size(z,2)
    jbest= 0;
    nbest= inf;
    outer= inf;
    
    % search nearest neighbour probabilistically
    for j=1:Nf
        [nis, nd]= compute_association(particle,z(:,i),R, j); % Nearest Neibour 
        if nis < gate1 && nd < nbest % if within gate, store nearest-neighbour
            nbest= nd;
            jbest= j;
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

function [nis, nd]= compute_association(particle,z,R,idf)
% return normalised innovation squared (ie, Mahalanobis distance) and normalised distance
[z_predict,~,~,Sf,~]= compute_jacobians(particle, idf, R);
v= z-z_predict;   % innovation
v(2)= pi_to_pi(v(2));

nis= v'/Sf*v;
nd= nis + log(det(Sf));



