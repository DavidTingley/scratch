% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');
binSize = [.01];
tcFactors = [2:2:20];
count_ls = 1;

for ii=60:length(d)
    cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    spikes = bz_GetSpikes('noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    if exist([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
    load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
    ls_spikesNREM = ls_spikes;
    ls_spikesBEHAV = ls_spikes;
    idx_ls = find(strcmp(spikes.region,'ls'));
    
    hpcRegion{ii} = 'ca1';
    if isempty(hpc_spikes) % comment out to exclude CA3 recs
        hpc_spikes = bz_GetSpikes('region','ca3','noprompts',true);
        hpcRegion{ii} = 'ca3';
    end
    idx_hpc = find(strcmp(spikes.region,'hpc') | strcmp(spikes.region,'ca3'));
    
    if exist([sessionInfo.FileName '.behavior.mat']) & ...
        ~isempty(ripples) & ...
        isfield(SleepState.ints,'NREMstate') & ~isempty(SleepState.ints.NREMstate)
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
    
    
    if ~isempty(ls_spikes)
    for spk = 1:length(ls_spikes.times)
%        ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},double(SleepState.ints.NREMstate)); 
       ls_spikesBEHAV.times{spk} = Restrict(ls_spikes.times{spk},behavior.events.trialIntervals);
       ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
    end
    end
   
    for bins = 1:length(binSize)
    if ~isempty(ls_spikes) & ~isempty(SleepState.ints.NREMstate) & exist([ls_spikes.sessionName '.behavior.mat']) & all(diff(intervals')>600)
        
    spkmat_ls = bz_SpktToSpkmat(ls_spikesBEHAV.times,'binSize',binSize(bins));
    spkmatNREM_ls = bz_SpktToSpkmat(ls_spikesNREM.times,'binSize',binSize(bins));
    
    for t = 1:length(firingMaps.rateMaps)
        template = squeeze(mean(firingMaps.rateMaps{t},2));
        figure(1)
        avgTemplateTime = median(diff(behavior.events.trialIntervals(behavior.events.trialConditions==t,:)'));
        [templateMatches pkRatio_low corrsAll] = bz_temporalCompression(spkmatNREM_ls,template(idx_ls,:)',avgTemplateTime,tcFactors);
        for tc = 1:length(tcFactors)
           preSleep_ls{count_ls}(tc,:) = corrsAll(tc,spkmatNREM_ls.timestamps<intervals(1,2));
           postSleep_ls{count_ls}(tc,:) = corrsAll(tc,spkmatNREM_ls.timestamps>intervals(3,1));
           behav_ls{count_ls}(tc,:) = corrsAll(tc,spkmatNREM_ls.timestamps>intervals(1,2) &...
                                spkmatNREM_ls.timestamps<intervals(3,1));
        end
        figure(2)
        subplot(3,2,1)
        plot(count_ls,nanmean(preSleep_ls{count_ls}(5,:)),'.g')
        hold on
        plot(count_ls,nanmean(postSleep_ls{count_ls}(5,:)),'.r')
        plot(count_ls,nanmean(behav_ls{count_ls}(5,:)),'.k')
        subplot(3,2,2)
        plot(count_ls,skewness(preSleep_ls{count_ls}(5,:)),'.g')
        hold on
        plot(count_ls,skewness(postSleep_ls{count_ls}(5,:)),'.r')
        plot(count_ls,skewness(behav_ls{count_ls}(5,:)),'.k')
        diff_ls(count_ls) = nanmean(preSleep_ls{count_ls}(5,:)) - nanmean(postSleep_ls{count_ls}(5,:));
        subplot(3,2,5)
        histogram(diff_ls)
        pause(.1)
        count_ls = 1+count_ls;
 
    end    
    end
    
    
    
    if ~isempty(hpc_spikes) 
    hpc_spikesNREM = hpc_spikes;
    hpc_spikesBEHAV = hpc_spikes;
    for spk = 1:length(hpc_spikes.times)
%         if length(hpc_spikes.times{spk})./hpc_spikes.times{spk}(end) < 5 % 2.5Hz FR limit
%        hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},double(SleepState.ints.NREMstate));   
       hpc_spikesBEHAV.times{spk} = Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals); 
       hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
%         else
%        hpc_spikesNREM.times{spk} = [];   
%        hpc_spikesBEHAV.times{spk} = []; 
%        disp('excluded a cell...')
%         end
    end
    spkmat_hpc = bz_SpktToSpkmat(hpc_spikesBEHAV.times,'binSize',binSize(bins));
    spkmatNREM_hpc = bz_SpktToSpkmat(hpc_spikesNREM.times,'binSize',binSize(bins));
    for t = 1:length(firingMaps.rateMaps)
        template = squeeze(mean(firingMaps.rateMaps{t},2));
        avgTemplateTime = median(diff(behavior.events.trialIntervals(behavior.events.trialConditions==t,:)'));
        [templateMatches pkRatio_low corrsAll] = bz_temporalCompression(spkmatNREM_hpc,template(idx_hpc,:)',avgTemplateTime,tcFactors);
        for tc = 1:length(tcFactors)
           preSleep_hpc{count_hpc}(tc,:) = corrsAll(tc,spkmatNREM_hpc.timestamps<intervals(1,2));
           postSleep_hpc{count_hpc}(tc,:) = corrsAll(tc,spkmatNREM_hpc.timestamps>intervals(3,1));
           behav_hpc{count_hpc}(tc,:) = corrsAll(tc,spkmatNREM_hpc.timestamps>intervals(1,2) &...
                                spkmatNREM_ls.timestamps<intervals(3,1));
            figure(2)
            subplot(3,2,3)
            plot(count_hpc,nanmean(preSleep_hpc{count_hpc}(5,:)),'.g')
            hold on
            plot(count_hpc,nanmean(postSleep_hpc{count_hpc}(5,:)),'.r')
            plot(count_hpc,nanmean(behav_hpc{count_hpc}(5,:)),'.k')
            subplot(3,2,4)
            plot(count_hpc,skewness(preSleep_hpc{count_hpc}(5,:)),'.g')
            hold on
            plot(count_hpc,skewness(postSleep_hpc{count_hpc}(5,:)),'.r')
            plot(count_hpc,skewness(behav_hpc{count_hpc}(5,:)),'.k')
            diff_hpc(count_hpc) = nanmean(preSleep_hpc{count_hpc}(5,:)) - nanmean(postSleep_hpc{count_hpc}(5,:));
            subplot(3,2,6)
            histogram(diff_hpc)
            pause(.1)
            count_hpc = 1+count_hpc;
        end       
    end
    end
    behavType{ii} = behavior.description;
   
    end
    end
    end
cd ~/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
end
