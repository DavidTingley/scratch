


d  = dir('*201*');
for ii=33:length(d)
   cd(d(ii).name) 
   spikes = bz_GetSpikes;
   load([spikes.sessionName '.behavior.mat'])
   
   load([spikes.sessionName '.sessionInfo.mat'])
   lfp = bz_GetLFP(sessionInfo.thetaChans(2));
   [rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,4);
   
    figure(1)
    clf
    bz_plotTrials(behavior)

    figure(2)
    f = factor(length(unique(behavior.events.trialConditions)));
    for i = 1:length(spikes.times)
        if strcmp(spikes.region{i},'ls')
            for j=1:length(unique(behavior.events.trialConditions))
                if ~isempty(phaseMap{j})
                    if ~isempty(phaseMap{j}{i})
            subplot(f(1),length(unique(behavior.events.trialConditions))./f(1),j)
            plot(1)
            scatter(phaseMap{j}{i}(:,1),phaseMap{j}{i}(:,end),'.k');
            hold on
            scatter(phaseMap{j}{i}(:,1),phaseMap{j}{i}(:,end)+2*pi,'.k');
            hold off
            axis([0 200 -pi pi*3])
                    end
                end
            end
            title(i)
            pause
        end
    end

    cd /home/david/datasets/lsDataset/
end