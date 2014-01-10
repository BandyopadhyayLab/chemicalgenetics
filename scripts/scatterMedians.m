%takes in the out object and scatter wt vs mut medians and labels points
%you can plot the data then make them in the color/opacity you want in
%illustrator. 
%has the ability to highlight certain spots of you tell it to
% raw_or_avg should be 1 if you want to work on the averaged data
function scatterMedians(out, raw_or_avg)


 if raw_or_avg==1
     out.scores = out.avgscores;
     out.names = out.avgnames;
     out.mag = out.avgmag;
     out.wtmed = out.avgwtmed;
     out.mutmed = out.avgmutmed;
 end
 figure
 hold
% scatter(out.wtmed(:),out.mutmed(:),'.')

 
 %scatter( out.wtmed(:) +
 %out.mutmed(:),log(out.mutmed(:)./out.wtmed(:)),'.')
 
 %scatter( out.wtmed(:) + out.mutmed(:),out.scores(:),'.')
 
 x=find(abs(out.scores)>=0 & abs(out.scores)<2);
 scatter(out.wtmed(x),out.mutmed(x),'b.')
 
  x=find(abs(out.scores)>=2 & abs(out.scores)<3);
 scatter(out.wtmed(x),out.mutmed(x),'y.')
 
  x=find(abs(out.scores)>=3 & abs(out.scores)<4);
 scatter(out.wtmed(x),out.mutmed(x),'g.')
   
 x=find(abs(out.scores)>=4 );
 scatter(out.wtmed(x),out.mutmed(x),'r.')
 
 legend('|S|<2','|S|=2-3','|S|=3-4','|S|>4')
 xlabel('Normalized cell number, control')
 
 ylabel('Normalized cell number, mutant')
 
 [i,j] = find(abs(out.scores)>=0);
 
 %for the negative interactions
 for index = [1:size(i,1)]
    
     
     loc = out.wtmed(i(index),j(index));
     loc2 =out.mutmed(i(index),j(index));
     
     name = out.names(i(index));
     label = out.druglabels(j(index));
     
     if size(cell2mat(strfind(name,'MYC')),1) && size(cell2mat(strfind(label,'Dasat')),1)
         
        text(loc,loc2,strcat(name,label),'FontSize',8)
        scatter(loc,loc2,'k')
     end
     
      if size(cell2mat(strfind(name,'MYC')),1) && size(cell2mat(strfind(label,'CHIR-99021')),1)
         
        text(loc,loc2,strcat(name,label),'FontSize',8)  
        scatter(loc,loc2,'k')
      end
       if size(cell2mat(strfind(name,'RAS')),1) && size(cell2mat(strfind(label,'Erlotinib')),1)
         
        text(loc,loc2,strcat(name,label),'FontSize',8) 
        scatter(loc,loc2,'k')
       end

 end

 
 figure
 hold

 axis2 = abs(out.scores);
 axis1 = out.mag;
 
 scatter( axis1(:),axis2(:),'.')
 
  %for the positive interactions

  for index = [1:size(i,1)]
    
     
     loc = axis1(i(index),j(index));
     loc2 =axis2(i(index),j(index));
     
     name = out.names(i(index));
     label = out.druglabels(j(index));

     if size(cell2mat(strfind(name,'MYC')),1) && size(cell2mat(strfind(label,'Dasat')),1)
         
        text(loc,loc2,strcat(name,label),'FontSize',8)
        scatter(loc,loc2,'k')
     end
     
      if size(cell2mat(strfind(name,'MYC')),1) && size(cell2mat(strfind(label,'CHIR-99021')),1)
         
        text(loc,loc2,strcat(name,label),'FontSize',8)  
        scatter(loc,loc2,'k')
      end
       if size(cell2mat(strfind(name,'RAS')),1) && size(cell2mat(strfind(label,'Erlotinib')),1)
         
        text(loc,loc2,strcat(name,label),'FontSize',8) 
        scatter(loc,loc2,'k')
       end
  end
  xlabel('Magnitude')
  ylabel('Score (|S|)')
  
end
