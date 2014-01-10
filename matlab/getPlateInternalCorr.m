%takes a 384 well plate and returns the average of the internal replicate
%correlations
function c=getPlateInternalCorr(plate)

if size(plate,1) == 0
    c=0;
    return 
end

if size(plate,1) ==16
    
    a=extract96from384(plate);
else
    a = plate;
end
tot = myNanCorrcoef(a(:,1),a(:,2))+ myNanCorrcoef(a(:,1),a(:,3))+ myNanCorrcoef(a(:,1),a(:,4))+myNanCorrcoef(a(:,2),a(:,3))+myNanCorrcoef(a(:,2),a(:,4))+myNanCorrcoef(a(:,3),a(:,4));


c=tot/6;


end
