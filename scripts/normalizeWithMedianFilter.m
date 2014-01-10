% Function to normalize out growth rate differences between control and
% oncogene expressing plates. Uses a sliding window to calculate the
% expected average cell number per well based on the assumption that most
% drugs to not have a differential effect on control versus oncogene
% expressing cell lines.  
%   wt = plate 1
%   mut = plate 2
%   plotson = 1 if you want to display the normalization details
%   numbins = value for the number of bins to use in this analysis.
% Author: SB
function [newwt,newmut]=normalizeWithMedianFilter(wt,mut,plotson,numbins)

    
    bins = numbins;
    %defines the limits of each evenly spaced bin based on min and max cell
    %number
    ind = [min(wt(:)):(max(wt(:))-min(wt(:)))/bins :max(wt(:))];
    
    wtmedians = ones(bins,1);
    mutmedians = ones(bins,1);

    newmut = mut;
    
    %special case for the first bin
    i = 1;
    %get the median of all values within that bin and two adjacent bins
    wtmedians(i) = nanmedian(wt(find(wt >= ind(i) & wt <= ind(i+2))));
    %find those same points on the mutant plate and adjust their median to
    %match with the control plate by adding an appropriate value.
    vals = mut(find(wt >= ind(i) & wt <= ind(i+2)));
    indexes  = find(wt >= ind(i) & wt <= ind(i+2));
    mutmedians(i)= nanmedian(vals);
    off=wtmedians(i)-mutmedians(i);
    newmut(indexes) = mut(indexes)+off;

       
    %all following bins
    for i=[2:bins-1]

        wtmedians(i) = nanmedian(wt(find(wt >= ind(i-1) & wt <= ind(i+1))));
        vals = mut(find(wt >= ind(i-1) & wt <= ind(i+1)));
        indexes  = find(wt >= ind(i-1) & wt <= ind(i+1));
        mutmedians(i)= nanmedian(vals);
        off=wtmedians(i)-mutmedians(i);
        newmut(indexes) = mut(indexes)+off;

        
    end

    %special case for the last bin
    vals = mut(find(wt>=ind(bins-2) & wt<=max(wt(:))));
    mutmedians(bins) = nanmedian(vals);
    wtmedians(bins) = nanmedian(wt(find(wt >= ind(bins-2) & wt <= max(wt(:)))));
    indexes  = find(wt >= ind(bins-2) & wt <= max(wt(:)));
    off=wtmedians(bins)-mutmedians(bins);
    newmut(indexes) = mut(indexes)+off;
    

    newwt = wt;
    if plotson == 1
        plot(wt,mut,'ro')
        hold
        plot(newwt,newmut,'bo')

        figure
        hold
        plot(newwt,newmut,'bo')
    end
end
