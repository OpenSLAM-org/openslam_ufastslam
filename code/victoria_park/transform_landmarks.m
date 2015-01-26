f=[];
k=0;
temp=[];
a=[];

for i=1:size(landmarks,2) % time step
    if landmarks(1,i) ~= landmarks(2,i)
        temp= [temp landmarks(:,i)]; k=k+1;
    else
        a=[a; landmarks(1,i)];
        
        %% Size adjusting
        d=size(f,2) - size(temp,2);
        
        
        if d > 0
            if d == size(f,2)
                temp=zeros(3,d);
            else
                temp=[temp zeros(3,d)];  
            end
        %%%
        elseif d < 0
            f=[f zeros(size(f,1),abs(d))];
        end
        %%%%%%%
        
        f=[f;temp];        
        
        k=0; temp=[];
    end
end

d=size(f,2) - size(temp,2);
if d > 0
    temp=[temp zeros(3,d)];
elseif d < 0
    f=[f zeros(size(f,1),abs(d))];
end
%%%%%%%
f=[f;temp];        
