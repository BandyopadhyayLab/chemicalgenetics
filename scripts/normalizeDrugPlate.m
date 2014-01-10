%converts the plate into a replicate noramlized z-scores
% deleted outliers at the 95% level using  the grubbs test method
%returns both the 96x4 and 384well shape versions

function shape96=normalizeDrugPlate(plate)
    


%remove these NaN replacements if you are using a custom drug library
% they were placed to remove drugs which were autoflourescent in Martins et
% al 2014.

    plate(5,9) = NaN;
    plate(5,10) = NaN;
    plate(6,9) = NaN;
    plate(6,10) = NaN;
    plate(11,23) = NaN;
    plate(11,24) = NaN;
    plate(12,23) = NaN;
    plate(12,24) = NaN;

    plate(13,19) = NaN;
    plate(13,20) = NaN;
    plate(14,19) = NaN;
    plate(14,20) = NaN;
    
    plate(7,15) = NaN;
    plate(7,16) = NaN;
    plate(8,15) = NaN;
    plate(8,16) = NaN;
    
    
   
    
    

    d1=extract96from384(plate);    
    
    norm_d1 = d1;
    
    for i =[1:96]
      
        if (min(isnan(norm_d1(i,:))) <1)
            
            norm_d1(i,:) = deleteoutliers(norm_d1(i,:),0.01,1);
        end
    end
    shape96 = norm_d1;
    
    
end


