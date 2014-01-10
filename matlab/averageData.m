%averages replicate scores and filters data where the variance between
%replicates is > varcutoff
function newout = averageData(out, varcutoff)
    newout = out;
    %get all the unique names
    for i = [1:size(out.names,1)]
    
        x=char(strsplit(char(out.names(i)),'_'));
        names(i) = cellstr(x(1,:));
    end

     
     unames = unique(names)';
     
     newout.avgnames = unames;
          
     for i = [1:size(unames,1)]
       ind= strmatch(unames(i),names);
     
       
       
       newout.avgscores(i,:) = nanmean(out.scores(ind,:));
       
       newout.avgmag(i,:) = nanmean(out.mag(ind,:));
       newout.avgwtvar(i,:) = nanmean(out.wtvar(ind,:));
       newout.avgmutvar(i,:) = nanmean(out.mutvar(ind,:));
       newout.avgwtmed(i,:) = nanmean(out.wtmed(ind,:));
       newout.avgmutmed(i,:) = nanmean(out.mutmed(ind,:));
        
       
       %remove scores that had high variance
       % we use 4
       newout.varscores(i,:) = nanvar(out.scores(ind,:));
       to_make_nan = find(newout.varscores(i,:)>=varcutoff);
       
       newout.avgscores(i,to_make_nan) = nan;
       newout.avgmag(i,to_make_nan) = nan;
       
     end
     
     
     
end
