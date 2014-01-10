%this function takes all the plates with the same gene and scatters their
%scores against each other and outputs their correlation
%also plots a confidence curve
% out = output from the removePlatesWithMinimumInternalCorr function
% Author: SB
function r=scatterReplicatePlates(out)
    c_array = [];
    r_array = [];
    local_r = [];
    badcorrcount = zeros(1,size(out.names,1));
    for i = [1:size(out.names,1)]
        if ~ismember(i,out.removedplates)
            a=out.names(i);
            name_one = regexp(a{1},'\_','split');

            for j=[i+1:size(out.names,1)]
                
                if ~ismember(j,out.removedplates)
                    a=out.names(j);
                    name_two = regexp(a{1},'\_','split');

                    %matching plates
                    if strcmp(name_one(1),name_two(1)) == 1
                       
                        
                       o = [out.scores(i,:)' , out.scores(j,:)']; 
                       
                       local_r(i) =  myNanCorrcoef(o(:,1),o(:,2));
                       %identification of plates that are not correlating
                       %with their biological replicates
                       if local_r(i)<0
                        
                        name_one;
                        name_two;
                        local_r(i);
                        
                        badcorrcount(i) = badcorrcount(i)+1;
                        badcorrcount(j) =badcorrcount(i)+1;
                        
                       end
                           
                           
                       c_array = [c_array ; o];
                       
                    %not matching plates, non replicates
                    else
                       o = [out.scores(i,:)' , out.scores(j,:)']; 
                            
                       r_array = [r_array ; o];
                      
                        
                    end

                 end



            end
        end
        
    end
    
    figure  
    scatter(c_array(:,1),c_array(:,2),'.')
    
    %figure
    %scatter(r_array(:,1),r_array(:,2),'.')
    
    
    replicatecorrelation = myNanCorrcoef(c_array(:,1),c_array(:,2))
    %hist(local_r,50);
    title(replicatecorrelation)
    xlabel('Score replicate 1')
    ylabel('Score replicate 2')
    %list of plates that are negatively correlated with their biological
    %replicates
    badplates = out.names(find(badcorrcount>=2))
    
   
    
end
