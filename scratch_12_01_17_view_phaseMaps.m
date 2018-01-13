d = dir('*201*')
for rec = 51:length(d)
    cd(d(rec).name)
    load([d(rec).name '.phaseMaps.cellinfo.mat'])
    load([d(rec).name '.behavior.mat'])
    
    figure(1)
    clf
    bz_plotTrials(behavior)
    
    figure(2)
    clf
    for cell=1:length(phaseMaps.UID)
        if strcmp(phaseMaps.region{cell},'ls')
    for cond = 1:length(phaseMaps.phaseMaps)
        f = factor(length(unique(behavior.events.trialConditions)));
        subplot(f(1),length(unique(behavior.events.trialConditions))./f(1),cond)
        if ~isempty(phaseMaps.phaseMaps{cond}{cell})
            scatter(phaseMaps.phaseMaps{cond}{cell}(:,1),phaseMaps.phaseMaps{cond}{cell}(:,end),'.k');
            hold on
            scatter(phaseMaps.phaseMaps{cond}{cell}(:,1),phaseMaps.phaseMaps{cond}{cell}(:,end)+2*pi,'.k');
            axis([0 200 -pi pi*3])
            hold off
        end
    end
    pause
        end
    end
    cd ..
end