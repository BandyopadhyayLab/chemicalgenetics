% Function to calculate S-scores from control and mutant drug response plates.
%   filename = name of the gene (filename) to analyze
%   wt = control plate data
%   mut = mutant cell line data
%   var_wt = variance of each of 96 wells across all control plates
%   plotson = 1 to see the plots, 0 otherwise
%   annotation_dir = location of the "labels" file to plot the drug names.
% Author: SB
function [s_score, mag, mut_med, wt_med, mut_var, wt_var] = compareDrugPlates(filename,wt,mut, var_wt, plotson, annotation_dir)

%adjust median
d = nanmedian(wt(:)) - 2000;
wt = wt  - d;

d1m = extract96from384(wt);
d2m = extract96from384(mut);

%normalize using a sliding window, in both directions and average
[d1m_1,d2m_1] = normalizeWithMedianFilter(d1m,d2m,0, 25);
[d1m_2,d2m_2] = normalizeWithMedianFilter(d2m,d1m,0, 25);

d1m = (d1m_1+d2m_2) / 2;
d2m = (d2m_1+d1m_2) / 2;

[order,ind]= sort(nanmedian(d2m') -nanmedian(d1m'));

wtm = nanmedian(d1m');

%scoring function using a minimum bound on variance so to avoid spuriously
%high signtificance values.
[h,p]=ttest_knownvariance(d1m',d2m',var_wt,nanvar(d2m'),'','','unequal','',20000);


m = nanmedian(d2m') - wtm;

st = nanstd(d2m'-[wtm' wtm' wtm' wtm']');


s_score = (-1*log10(p)) .* sign(m);
v1 = ones(1,101);
x  = [0:100]/100;



mag = log2(nanmedian(d2m')  ./ nanmedian(d1m'));
mut_med = nanmedian(d2m');
wt_med = nanmedian(d1m');
wt_var = nanvar(d1m');
mut_var = nanvar(d2m');

if (plotson == 1)
    labels =textread(strcat(annotation_dir,'labels'),'%s');

    figure
    hold
    scatter(nanmedian(d1m') ,  nanmedian(d2m'),'.')
    errorbar(nanmedian(d1m') ,  nanmedian(d2m'), nanstd(d2m')./2 , 'o');
    herrorbar(nanmedian(d1m') ,  nanmedian(d2m'), nanstd(d1m')./2, 'o');

    text(nanmedian(d1m') ,  nanmedian(d2m'),labels);

    i=find(p<0.001);
    scatter(nanmedian(d1m(i,:)') ,  nanmedian(d2m(i,:)') ,'r','filled');
    xlabel('Control');
    ylabel(filename);
    figure

    scatter(mag,abs(s_score),'.')


    i = union(find(mag>0.5 | mag<-0.5),find(abs(s_score)>1.5));
    
    sigmag = nanmedian(d2m(i,:)')  ./ nanmedian(d1m(i,:)') -1;
    % make the max mag 4 for plotting purposes
    sigmag(find(sigmag>=4)) = 4;
    sigmag(find(sigmag<=-4)) = -4;

    text(sigmag,abs(s_score(i)),labels(i))
end



