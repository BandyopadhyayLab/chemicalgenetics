%processes the MCF10A single drug screen by batch, compared to controls in
% directory
%   directory_name = the location of the experiment plates 
%   controls_dir = the location of the control plates
%   numbatches = the number of batches that the plates were collected in corresponding to the maximum B# in the file names.
%   annotation_dir = the location of the "labels" file containing the names
%       of the drugs.
%   plotson = 1 if you would like to see the graph based analyses, 0
%       otherwise
% Author: SB
function [output]=loadDirectoryofPlates(directory_name, controls_dir, numbatches, annotation_dir, plotson)

    files = dir(controls_dir);
    fileIndex = find([files.bytes]>0);
    files2 = dir(directory_name);
    fileIndex2 = find([files2.bytes]>0);
    n = length(fileIndex2);
    k=cell(n,1);
    output.names = k;
    output.scores = ones(length(fileIndex2),96)-1;
    labels =textread(strcat(annotation_dir,'labels'),'%s');
    
    % maximum CV that is allowed among the 4 replicates across 96 samples in the 384 well
    % plate.
    cv_cutoff = 0.4;
    
    for batch = [1:numbatches]
        st_control = strcat('_B',num2str(batch),'_');
        st = strcat('_B',num2str(batch),'_');
        

        mergedplate = ones(16,24)-1;
        numels = ones(16,24)-1;
        plate_aggregate = [];
        for i = 1:length(fileIndex)
                
                fileName = files(fileIndex(i)).name;
                if strfind(fileName,st_control)>1
                    fileName
                    current = load(strcat(controls_dir,fileName));
                    
                    %remove outliers and remove problematic drugs
                    out = normalizeDrugPlate(current);
                    
                    current = extract384from96(out);
                    %adjust median to 2000
                    current = current * (2000/nanmedian(current(:)));
                    plate_aggregate = [plate_aggregate extract96from384(current)];
                    for index = [1:384]
                        if isnan(current(index))==0
                            mergedplate(index) = mergedplate(index)+current(index);
                            numels(index) = numels(index)+1;
                        end

                    end
                    
                    
                    
                end
                
        end
 
                
        %go through all the control plates in this batch and
        %determine the CV 
        cv = nanstd(plate_aggregate') ./ nanmean(plate_aggregate');
        %for drugs that have CV higher than cutoff remove the drug from
        %further analysis
         if size(cv,2)>1
             temp = extract384from96([cv' cv' cv' cv']);
             toremove = find(temp>=cv_cutoff);
 
             mergedplate(toremove) = nan;
         end
        
        mergedplate = mergedplate ./ numels;


        num_drugs_removed = size(find(cv>=cv_cutoff))
        
        
        for queryplate = 1:length(fileIndex2)
            fileName = files2(fileIndex2(queryplate)).name;
            if strfind(fileName,st)>1
                   
                current = load(strcat(directory_name,fileName));
                %remove outliers
                out = normalizeDrugPlate(current);
                current = extract384from96(out);
                %adjust median to 2000
                current = current * (2000/nanmedian(current(:)));

                fileName
                [s_score,mag,mutmed, wtmed, mutvar, wtvar]= compareDrugPlates(fileName,mergedplate,current,nanvar(plate_aggregate'),plotson, annotation_dir);
                
                %get the normalized plate count data 
                d1m = normalizeDrugPlate(mergedplate);
                d2m = normalizeDrugPlate(current);
                [d1m,d2m] = normalizeWithMedianFilter(mergedplate,current,0, 5);
                
                ic = getPlateInternalCorr(current);
                %do not analyze plates with this minimum internal
                %correlations
                if ic>0.7
                    output.names(queryplate) = java.lang.String(fileName);
                    output.scores(queryplate,:) = s_score();
                    output.mutvar(queryplate,:) = mutvar();
                    output.wtvar(queryplate,:) = wtvar();
                    output.mutmed(queryplate,:) = mutmed();
                    output.wtmed(queryplate,:) = wtmed();
                    output.queryarray(queryplate,:) = d2m(:);
                    output.controlarray(queryplate,:) = d1m(:);
                    output.internalcorr(queryplate,:) = ic;
                    output.mag(queryplate,:) = mag();

                end

            close all;
        end
    end
    
   end
    
    output.removedplates = [];
    output.druglabels = labels();
    %amount of data points that were removed because of 
    percent_missing_data = sum(isnan(output.scores(:))) ./ size(output.scores(:),1)
      
    
    % removes any plates where the correlation of scores between replicates
    % is <0.7
    out2 = removePlatesWithMinimumInternalCorr(output,0.7); 
    
    %scatter replicate genetic interaction scores.
    scatterReplicatePlates(out2);
    out3 = averageData(out2,4);
    scatterMedians(out3,1)    
    %return data
    output = out2;
     
    
    
     
    
    