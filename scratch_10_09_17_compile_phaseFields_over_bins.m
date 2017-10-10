















clear all
d = dir('*201*');
count=1;
for rec = 1:length(d)
    cd(d(rec).name)
    load([d(rec).name '.firingMaps.cellinfo.mat'])
    load([d(rec).name '.behavior.mat'])
    if strcmp(behavior.description,'central alternation')
    for cell=1:length(firingMaps.UID)
    if strcmp(firingMaps.region{cell},'ls')
       for cond = 1:length(firingMaps.rateMaps) 
           if ~isempty(firingMaps.phaseMaps{cond}{cell})
        for bin = 1:200
            f = find(abs(firingMaps.phaseMaps{cond}{cell}(:,1)-bin) < 20);
            if length(f) > 20
                [pval(count) z] = circ_rtest(firingMaps.phaseMaps{cond}{cell}(f,1));
                [rval(count)] = circ_r(firingMaps.phaseMaps{cond}{cell}(f,1));
                [mval(count)] = circ_mean(firingMaps.phaseMaps{cond}{cell}(f,1));
                [sval(count)] = circ_std(firingMaps.phaseMaps{cond}{cell}(f,1));
                bins(count) = bin;
                count=1+count;
            end
        end
           end 
       end
    end
    end
    end
    cd /home/david/datasets/lsDataset
end

clf
for i=1:200
f = find(bins==i & pval < .05);
scatter(i,mean(rval(f)'),'.k')
hold on
end