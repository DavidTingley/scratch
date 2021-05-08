clear all
orig_files = dir('_*mat');
new_files = dir('revision/*mat');

cgm_files = dir('stim/cgm*mat');
% cgm_files = [orig_files;new_files];
troughTrig_all = [];
peakTrig_all = [];

for j=1:length(cgm_files)
    if ~strcmp(cgm_files(j).name,'CGM36.mat')
        load([cgm_files(j).folder '/' cgm_files(j).name],'count','stimRate','isig_levels','isa','idx','spSlope','emgSig','mov','theta_z','states')
        count = stimRate; 
        if ~strcmp(cgm_files(j).name,'CGM47.mat') & ~strcmp(cgm_files(j).name,'CGM50.mat')
            isig = nanZscore([0, diff(isig_levels)]);
        else
            isig = nanZscore(isa);
        end
        if strcmp(cgm_files(j).name,'ros.mat')
            isig = nanZscore(isa);
        end

        [pks locs] = findpeaks(-isig,'MinPeakHeight',1.5);
        locs = locs(find(ismember(locs,idx)));
        for i=1:length(locs)
            if locs(i)>40 & locs(i) < length(count)-40
            rr(i,:)=count(locs(i)-40:locs(i)+40);
            end
        end
        
        
        [pks locs] = findpeaks(isig,'MinPeakHeight',1.1);
        locs = locs(find(ismember(locs,idx)));
        for i=1:length(locs)
            if locs(i)>40 & locs(i) < length(count)-40
            r(i,:)=count(locs(i)-40:locs(i)+40);
            end
        end
        
        if exist('r')
            peakTrig_all = [peakTrig_all;r];
            peakTrig(j,:) = nanmean(zscore(r,[],2));
            for it = 1:100
            peakTrig_s(j,it,:) = nanmean(zscore(bz_shuffleCircular(r),[],2)); 
            end
        end
        clear r;
        if exist('rr')
            troughTrig_all = [troughTrig_all;rr];
            troughTrig(j,:) = nanmean(zscore(rr,[],2)); 
            for it = 1:100
            troughTrig_s(j,it,:) = nanmean(zscore(bz_shuffleCircular(rr),[],2)); 
            end
        end
        clear rr;
    end
end

clf
boundedline(1:81,squeeze(nanmean(nanmean(troughTrig_s))),squeeze(nanstd(nanmean(troughTrig_s)))*3,'r')
boundedline(1:81,squeeze(nanmean(nanmean(peakTrig_s))),squeeze(nanstd(nanmean(peakTrig_s)))*3,'k')
xline(41)
hold on
plot(nanmean(peakTrig),'k')
plot(nanmean(troughTrig),'r')
title('red/trough; black/peaks')
ylabel('z-scored rate')