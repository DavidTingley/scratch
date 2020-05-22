cd D:\datasets\ripples_LS
d = dir('*201*');

for i =68:98
    cd(d(i).name)
    if exist([d(i).name '.phaseMaps.cellinfo.mat'])
    load([d(i).name '.phaseMaps.cellinfo.mat']);
    load([d(i).name '.behavior.mat']);
    spikes = bz_GetSpikes('noprompts',true);
    
    
    for cond = 1:length(behavior.events.map)
        for s=1:length(spikes.times)
        spk.times{s} = Restrict(spikes.times{s},[behavior.events.trialIntervals(behavior.events.trialConditions==cond,:)]);
        spk_all.times{s} = Restrict(spikes.times{s},[behavior.events.trialIntervals]);
        end
        [times groups] = spikes2sorted(spk_all.times);
        [ccg_beh_all t] = CCG(times,groups,'binSize',.001,'duration',.5);
        
        [times groups] = spikes2sorted(spk.times);
        [ccg_beh t] = CCG(times,groups,'binSize',.001,'duration',.5);
    
        
        for s=1:length(spikes.times)
            for ss = 1:length(spikes.times)
                if ss>s & strcmp(spikes.region{s},'ls') & strcmp(spikes.region{ss},'ls') & sum(ccg_beh(:,s,ss)) > 10
                    subplot(3,2,1)
                    scatter(phaseMaps.phaseMaps{cond}{s}(:,1),phaseMaps.phaseMaps{cond}{s}(:,end),'.k')
                    hold on
                    scatter(phaseMaps.phaseMaps{cond}{s}(:,1),phaseMaps.phaseMaps{cond}{s}(:,end)+2*pi,'.k')
                    hold off
                    title(s)
                    
                    subplot(3,2,2)
                    scatter(phaseMaps.phaseMaps{cond}{ss}(:,1),phaseMaps.phaseMaps{cond}{ss}(:,end),'.k')
                    hold on
                    scatter(phaseMaps.phaseMaps{cond}{ss}(:,1),phaseMaps.phaseMaps{cond}{ss}(:,end)+2*pi,'.k')
                    hold off
                    title(ss)
                    
                    subplot(3,2,3)
                    plot(ccg_beh_all(:,s,ss),'k')
                    hold on
                    plot(ccg_beh(:,s,ss),'r')
                    hold off
                    
                    subplot(3,2,4)
                    plot(Smooth(ccg_beh_all(:,s,ss),5),'k')
                    hold on
                    plot(Smooth(ccg_beh(:,s,ss),5),'r')
                    hold off
                    
                    subplot(3,2,5)
                    plot(ccg_post{i}(:,s,ss),'b')
                    hold on
                    plot(ccg_pre{i}(:,s,ss),'g')
                    hold off
                    title([excess_pre{i}(s,ss) excess_post{i}(s,ss)])
                    
                    pause
                end
            end
        end
    end
    
    end
    cd ..
end