% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');

for int =1:3
    ccg_ls{int} = [];
    ccg_hpc{int} = [];
    ccg_cross{int} = [];
end


for ii=61:length(d)
    cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    ls_spikesNREM = ls_spikes;
    ls_spikesBEHAV = ls_spikes;
    
    
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
    if ~isempty(ls_spikes) & ~isempty(SleepState.ints.NREMstate) & exist([ls_spikes.sessionName '.behavior.mat']) & all(diff(intervals')>600)
      
    for spk = 1:length(ls_spikes.times) 
       ls_spikesBEHAV.times{spk} = [0;Restrict(ls_spikes.times{spk},behavior.events.trialIntervals)];
       ls_spikesNREM_pre.times{spk} = [0;Restrict(Restrict(Restrict(ls_spikes.times{spk},[ripples.peaks-.2 ...
           ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
       ls_spikesNREM_post.times{spk} = [0;Restrict(Restrict(Restrict(ls_spikes.times{spk},[ripples.peaks-.2 ...
           ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
    end
    
%      for spk = 1:length(ls_spikes.times)
%        ls_spikesBEHAV.times{spk} = [0;Restrict(ls_spikes.times{spk},behavior.events.trialIntervals)];
%        ls_spikesNREM_pre.times{spk} = [0;Restrict(Restrict((ls_spikes.times{spk} ...
%            ),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
%        ls_spikesNREM_post.times{spk} = [0;Restrict(Restrict((ls_spikes.times{spk} ...
%            ),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
%     end
    end
   


    for int = 1:size(intervals,1)
        if int == 1 
            [times groups]= spikes2sorted(ls_spikesNREM_pre.times);
        elseif int == 3
            [times groups]= spikes2sorted(ls_spikesNREM_post.times);
        elseif int == 2
            [times groups]= spikes2sorted(ls_spikesBEHAV.times);
        end
        [ccg t] = CCG(times,groups,'binSize',.001);
        for c1 = 1:size(ccg,2)
            for c2 = c1:size(ccg,3)
                ccg_ls{int} = [ccg_ls{int}, zscore(Smooth(squeeze(ccg(:,c1,c2)),5))];
            end
        end
    end 
    clear ls_spikes*
    
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    hpcRegion{ii} = 'ca1';
    if isempty(hpc_spikes) % comment out to exclude CA3 recs
        hpc_spikes = bz_GetSpikes('region','ca3','noprompts',true);
        hpcRegion{ii} = 'ca3';
    end
    if ~isempty(hpc_spikes) 
    for spk = 1:length(hpc_spikes.times)
       hpc_spikesBEHAV.times{spk} = [0;Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals)];
       hpc_spikesNREM_pre.times{spk} = [0;Restrict(Restrict(Restrict(hpc_spikes.times{spk},[ripples.peaks-.2 ...
           ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
       hpc_spikesNREM_post.times{spk} = [0;Restrict(Restrict(Restrict(hpc_spikes.times{spk},[ripples.peaks-.2 ...
           ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
    end
    % no ripple restriction
%     for spk = 1:length(hpc_spikes.times)
%        hpc_spikesBEHAV.times{spk} = [0;Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals)];
%        hpc_spikesNREM_pre.times{spk} = [0;Restrict(Restrict((hpc_spikes.times{spk} ...
%            ),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
%        hpc_spikesNREM_post.times{spk} = [0;Restrict(Restrict((hpc_spikes.times{spk} ...
%            ),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
%     end

    for int = 1:size(intervals,1)
        if int == 1 
            [times groups]= spikes2sorted(hpc_spikesNREM_pre.times);
        elseif int == 3
            [times groups]= spikes2sorted(hpc_spikesNREM_post.times);
        elseif int == 2
            [times groups]= spikes2sorted(hpc_spikesBEHAV.times);
        end
        [ccg t] = CCG(times,groups,'binSize',.001);
        for c1 = 1:size(ccg,2)
            for c2 = c1:size(ccg,3)
                ccg_hpc{int} = [ccg_hpc{int}, zscore(Smooth(squeeze(ccg(:,c1,c2)),5))];
            end
        end
    end
    clear hpc_spikes*
    end
    
    %% cross region
    spikes = bz_GetSpikes('noprompts',true);
    if ~isempty(spikes) 
    for spk = 1:length(spikes.times)
       spikesBEHAV.times{spk} = [0;Restrict(spikes.times{spk},behavior.events.trialIntervals)];
       spikesNREM_pre.times{spk} = [0;Restrict(Restrict(Restrict(spikes.times{spk},[ripples.peaks-.2 ...
           ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
       spikesNREM_post.times{spk} = [0;Restrict(Restrict(Restrict(spikes.times{spk},[ripples.peaks-.2 ...
           ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
    end
    % no ripple restriction
%     for spk = 1:length(spikes.times)
%        spikesBEHAV.times{spk} = [0;Restrict(spikes.times{spk},behavior.events.trialIntervals)];
%        spikesNREM_pre.times{spk} = [0;Restrict(Restrict((spikes.times{spk} ...
%            ),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
%        spikesNREM_post.times{spk} = [0;Restrict(Restrict((spikes.times{spk} ...
%            ),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
%     end

    for int = 1:size(intervals,1)
        if int == 1 
            [times groups]= spikes2sorted(spikesNREM_pre.times);
        elseif int == 3
            [times groups]= spikes2sorted(spikesNREM_post.times);
        elseif int == 2
            [times groups]= spikes2sorted(spikesBEHAV.times);
        end
        [ccg t] = CCG(times,groups,'binSize',.001);
        for c1 = 1:size(ccg,2)
            for c2 = c1:size(ccg,3)
                if strcmp(spikes.region{c1},'ls') 
                    if strcmp(spikes.region{c2},'hpc') | strcmp(spikes.region{c2},'ca3')
                ccg_cross{int} = [ccg_cross{int}, zscore(Smooth(squeeze(ccg(:,c1,c2)),5))];
                    end
            end
        end
    end
    
    end
    clear spikes*
    behavType{ii} = behavior.description;
      

    pause(.1)

    end
    end
    
    
for warp = 1:10:999
    for c=1:size(ccg_hpc{1},2)
    cc(:,c) = makeLength(ccg_hpc{1}(warp:2001-warp,c),2001);
    ccc(:,c) = makeLength(ccg_hpc{3}(warp:2001-warp,c),2001);
    end
    pre_hpc(warp) = corr2(ccg_hpc{2},cc);
    post_hpc(warp) = corr2(ccg_hpc{2},ccc);
    pre_hpc_shuf(warp) = corr2(ccg_hpc{2},bz_shuffleCircular(cc));
    post_hpc_shuf(warp) = corr2(ccg_hpc{2},bz_shuffleCircular(ccc));
    subplot(2,2,1)
    plot(pre_hpc,'.k')
    hold on
    plot(pre_hpc_shuf,'.r')
    hold off
    subplot(2,2,2)
    plot(post_hpc,'.k')
    hold on
    plot(post_hpc_shuf,'.r')
    hold off
    title('HPC CCG reactivation')
    clear cc ccc
    
    for c=1:size(ccg_ls{1},2)
    cc(:,c) = makeLength(ccg_ls{1}(warp:2001-warp,c),2001);
    ccc(:,c) = makeLength(ccg_ls{3}(warp:2001-warp,c),2001);
    end
    pre_ls(warp) = corr2(ccg_ls{2},cc);
    post_ls(warp) = corr2(ccg_ls{2},ccc);
    pre_ls_shuf(warp) = corr2(ccg_ls{2},bz_shuffleCircular(cc));
    post_ls_shuf(warp) = corr2(ccg_ls{2},bz_shuffleCircular(ccc));
    subplot(2,2,3)
    plot(pre_ls,'.m')
    hold on
    plot(pre_ls_shuf,'.r')
    hold off
    subplot(2,2,4)
    plot(post_ls,'.m')
    hold on
    plot(post_ls_shuf,'.r')
    hold off
    title('LS CCG reactivation')
    clear cc ccc
    
    pause(.01)
end

    end
cd ~/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
   
end
    
    
    
    
    
