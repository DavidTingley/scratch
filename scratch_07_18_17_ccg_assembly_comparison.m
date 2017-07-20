


xml = LoadParameters;
load([xml.FileName '.sessionInfo.mat'])
load([xml.FileName '.spikes.cellinfo.mat'])
load([xml.FileName '.behavior.mat'])
load('assembliesCrossRegionData.mat')
lfp = bz_GetLFP(sessionInfo.thetaChans(2));
[rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,4);

for cond = 1:length(unique(behavior.events.trialConditions))
    for cell = 1:length(spikes.times)
        sp{cell} = Restrict(spikes.times{cell},behavior.events.trialIntervals(behavior.events.trialConditions==cond,:));
    end
    [times groups] = spikes2sorted(sp);
    [ccg t]  = CCG(times,groups,'binSize',.001,'duration',3);
    if cond <= length(dev)
    for pair = 1:size(dev{cond},2)
    [a b] =  min(dev{cond}(:,pair));
    [aa bb] = min(mean(devControl{cond}(:,pair,:),3));
    imp = (a-mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(:,pair,:),3)));
    imp2 = a ./ max(mean(mean(devControl{cond}(:,pair,:),3)));
    zerolag = (min(dev{cond}(1:6,pair)) - mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(1,pair,:),3)));
    if zerolag < 1 
        zerolag = 1;
    end
    
    if imp > 5 & b > 7 & zerolag < 1.2 & mean(dev{cond}(:,pair))>150
        ls = pairs(pair,1);
        hpc = pairs(pair,2);
        subplot(3,2,1)
        imagesc(squeeze(rateMap{cond}(ls,:,:)))
        title(ls)
        subplot(3,2,2)
        imagesc(squeeze(rateMap{cond}(hpc,:,:)))
        title(hpc)
        subplot(3,2,3)
        plot(t,smooth(ccg(:,ls,hpc),20))
        subplot(3,2,4)
        plot(dev{cond}(:,pair))
        hold on
        plot(squeeze(mean(devControl{cond}(:,pair,:),3)),'r');
        title(imp)
        hold off
        subplot(3,2,5)
        plot(spikes.rawWaveform{ls}*.195)
        hold on
        plot(spikes.rawWaveform{hpc}*.195,'r')
        hold off
        
        pause
    end
    end
end
end