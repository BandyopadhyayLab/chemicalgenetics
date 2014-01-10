function out = removePlatesWithMinimumInternalCorr(data, cutoff)

    out = data;
    out.orig_names = out.names;
    toremove = find(out.internalcorr<cutoff);
    out.names(toremove)= '';
    out.removedplates = toremove;

    out.scores = out.scores(setdiff(1:size(out.scores,1),toremove),:);
    out.mutvar = out.mutvar(setdiff(1:size(out.mutvar,1),toremove),:);
    out.mutmed = out.mutmed(setdiff(1:size(out.mutmed,1),toremove),:);
    out.wtvar = out.wtvar(setdiff(1:size(out.wtvar,1),toremove),:);
    out.wtmed = out.wtmed(setdiff(1:size(out.wtmed,1),toremove),:);
    out.mag = out.mag(setdiff(1:size(out.mag,1),toremove),:);

end
