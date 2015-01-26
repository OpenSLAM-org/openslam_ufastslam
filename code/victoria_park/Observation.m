%% laserfeature --> proper dimension
clear;
load 'feature.txt';
ObservationDistance=[];
ObservationAngle=[];

for i=1:3:size(feature,1)
    ObservationDistance = [ObservationDistance ; feature(i,:)];
end
for i=2:3:size(feature,1)
    ObservationAngle = [ObservationAngle ; feature(i,:)];
end

save 'ObservationDistance.txt' -ascii ObservationDistance;
save 'ObservationAngle.txt' -ascii ObservationAngle;

