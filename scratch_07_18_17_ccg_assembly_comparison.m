d  = dir('*201*');

for i=17:length(d)
cd(d(i).name) 

if exist('assembliesCrossRegionData_w_theta_sin_cos_coord_vel.mat')
xml = LoadParameters;
load([xml.FileName '.sessionInfo.mat'])
load([xml.FileName '.spikes.cellinfo.mat'])
load([xml.FileName '.behavior.mat'])
load([xml.FileName '.placeFields.20_pctThresh.mat'])
load([xml.FileName '.firingMaps.cellinfo.mat'])
load([xml.FileName '.phaseMaps.cellinfo.mat'])
try
    load('assembliesCrossRegion_split_w_theta_08-Nov-2017.mat')
catch
    load('assembliesCrossRegionData_w_theta_sin_cos_coord_vel.mat')
end
% lfp = bz_GetLFP(sessionInfo.thetaChans(2));
% [firingMaps] = bz_firingMap1D(spikes,behavior,lfp,4);

for cond = 1:length(unique(behavior.events.trialConditions))
    for cell = 1:length(spikes.times)
        sp{cell} = Restrict(spikes.times{cell},behavior.events.trialIntervals(behavior.events.trialConditions==cond,:));
    end
    [times groups] = spikes2sorted(sp);
    [ccg t]  = CCG(times,groups,'binSize',.001,'duration',1);
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
    if ~isempty(fields{cond}{pairs(pair,2)})
    if imp > 4 & b > 7 & b < 150 &  zerolag < 1.2 & mean(dev{cond}(:,pair))>50
        ls = pairs(pair,1);
        hpc = pairs(pair,2);
        subplot(4,2,1)
        imagesc(squeeze(firingMaps.rateMaps{cond}(ls,:,:)))
        title(ls)
        subplot(4,2,2)
        imagesc(squeeze(firingMaps.rateMaps{cond}(hpc,:,:)))
        
        title(hpc)
        subplot(4,2,3)
        
        scatter(phaseMaps.phaseMaps{cond}{ls}(:,1),(phaseMaps.phaseMaps{cond}{ls}(:,end)),'.k');
        hold on
        scatter(phaseMaps.phaseMaps{cond}{ls}(:,1),(phaseMaps.phaseMaps{cond}{ls}(:,end))+2*pi,'.k');
        title(firingMaps.sessionName)
        axis([0 200 -pi pi*3])
        hold off
        subplot(4,2,4)
        
        scatter(phaseMaps.phaseMaps{cond}{hpc}(:,1),(phaseMaps.phaseMaps{cond}{hpc}(:,end)),'.k');
        hold on
        scatter(phaseMaps.phaseMaps{cond}{hpc}(:,1),(phaseMaps.phaseMaps{cond}{hpc}(:,end))+2*pi,'.k');
        hold off
        axis([0 200 -pi pi*3])
        subplot(4,2,5)
        plot(spikes.rawWaveform{ls}*.195)
        hold on
        plot(spikes.rawWaveform{hpc}*.195,'r')
        hold off
        title('b-ls, r-hpc')
        subplot(4,2,6)
        plot(dev{cond}(:,pair))
        hold on
        plot(squeeze(mean(devControl{cond}(:,pair,:),3)),'r');
        title(imp)
        hold off
        
        subplot(4,2,7)
        plot(t,smooth(ccg(:,ls,hpc),20))
        title(['cond: ' num2str(cond)])
        pause
    end
    end
    end
end
end
end

cd /home/david/datasets/lsDataset
end