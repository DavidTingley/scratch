clear 
count = 1;
d = dir('*201*');

for rec=1:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    spikes = bz_GetSpikes('noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    if exist([sessionInfo.FileName '.behavior.mat']) & ...
        isfield(SleepState.ints,'NREMstate') & ~isempty(spikes) & ~isempty(SleepState.ints.NREMstate)
        load([sessionInfo.FileName '.behavior.mat'])
        for i=1:length(behavior.events.trials)
            bStart(i) = behavior.events.trials{i}.timestamps(1);
            bStop(i) = behavior.events.trials{i}.timestamps(end);
        end
        pre = [0 min(bStart)];
        post = [max(bStop)  sessionInfo.recordingDuration];
        behav = [pre(end) post(1)]; clear bStop bStart
        %% get initial correlations
        intervals = [pre; behav; post];
        for spk = 1:length(spikes.times)
            for int = 1:3
                if int == 1 | int == 3
                    sp = Restrict(Restrict(spikes.times{spk},intervals(int,:)),double(SleepState.ints.NREMstate));
                    rate(count,int) = length(sp)./sum(diff(SubtractIntervals(intervals(int,:),double(SleepState.ints.NREMstate))'));
%                     sp = Restrict(spikes.times{spk},intervals(int,:));
%                     rate(count,int) = length(sp)./diff(intervals(int,:));
                elseif int == 2
                    sp = Restrict(spikes.times{spk},behavior.events.trialIntervals);% intervals(int,:));
                    rate(count,int) = length(sp)./sum(diff(SubtractIntervals(intervals(int,:),behavior.events.trialIntervals)'));
                    
%                     sp = Restrict(spikes.times{spk},intervals(int,:));
%                     rate(count,int) = length(sp)./diff(intervals(int,:));
                end
                reg{count} = spikes.region{spk};
            end
            rate_z(count,:) = zscore(rate(count,:));
            rate_meanNorm(count,:) = (rate(count,:)./mean(rate(count,:)));
            
        count=1+count;
        end
    end
    cd ~/datasets/ripples_LS/
    
    if exist('rate_z')
        hpc = find(strcmp(reg,'hpc'));
        LS = find(strcmp(reg,'ls'));
        ca3 = find(strcmp(reg,'ca3'));
        subplot(3,2,1)
        errorbar(1:3,nanmedian(rate_z(hpc,:)),sem(rate_z(hpc,:)),'.r')
        hold on
        errorbar(1.1:3.1,nanmedian(rate_z(ca3,:)),sem(rate_z(ca3,:)),'.b')
        errorbar(1.2:3.2,nanmedian(rate_z(LS,:)),sem(rate_z(LS,:)),'.g')
        axis([0 4 -1.5 1.5])
        ylabel('zscored FR')
        xlabel('pre-NREM  behavior  post-NREM')
        hold off
        subplot(3,2,5)
        errorbar(1:3,nanmean(rate(hpc,:)),sem(rate(hpc,:)),'.r')
        hold on
        errorbar(1.1:3.1,nanmean(rate(ca3,:)),sem(rate(ca3,:)),'.b')
        errorbar(1.2:3.2,nanmean(rate(LS,:)),sem(rate(LS,:)),'.g')
        hold off
        subplot(3,2,6)
        errorbar(1:3,nanmean(rate_meanNorm(hpc,:)),3*sem(rate_meanNorm(hpc,:)),'.r')
        hold on
        errorbar(1.1:3.1,nanmean(rate_meanNorm(ca3,:)),3*sem(rate_meanNorm(ca3,:)),'.b')
        errorbar(1.2:3.2,nanmean(rate_meanNorm(LS,:)),3*sem(rate_meanNorm(LS,:)),'.g')
%         axis([0 4 -1.5 1.5])
        ylabel('FR')
        xlabel('pre-NREM  behavior  post-NREM')
        hold off
        
        subplot(3,2,2);
        imagesc(rate_z(LS,:))
        caxis([-2 2])
        title('ls')
        subplot(3,2,3)
        imagesc(rate_z(ca3,:))
        caxis([-2 2])
        title('ca3')
        subplot(3,2,4)
        imagesc(rate_z(hpc,:))
        caxis([-2 2])
        title('ca1')
        pause(.01)
    end
end



