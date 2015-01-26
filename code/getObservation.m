function z = getObservation(j,ObservationD,ObservationA,perceplimit)

z=[];

for i=1:size(ObservationD,2),
    if(ObservationD(j,i) ~= 0 && ObservationD(j,i) < perceplimit)
        z=[z [ObservationD(j,i);ObservationA(j,i)-pi*0.5]];
    end
end