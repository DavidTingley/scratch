% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');
binSize = [.025 1];%[.01 .015 .025 .05 .1 1 10];

for ii=1:length(d)
    cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    ripples = bz_LoadEvents(pwd,'LSRipples');
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
    wind = 1;
    for window = 0:300:3000
    % get the first XXX minutes of sleep..
    temp = 1;
    epoch = 0;
    while sum(epoch) < 600 & temp < intervals(end)
        it = SubtractIntervals(double(SleepState.ints.NREMstate),[0 intervals(2,2)+window; intervals(2,2)+temp+window intervals(end)]);
        epoch = sum(diff(it'));
        temp = temp+1;
    end
    
    
    
    if ~isempty(ls_spikes) & ~isempty(it)
    for spk = 1:length(ls_spikes.times)
       ls_spikesNREM_post.times{spk} = Restrict(ls_spikes.times{spk},it); 
%        ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},double(SleepState.ints.NREMstate)); 
       ls_spikesBEHAV.times{spk} = Restrict(ls_spikes.times{spk},behavior.events.trialIntervals);
       ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
    end
    end
   
    for bins = 1:length(binSize)
    if ~isempty(ls_spikes) & ~isempty(SleepState.ints.NREMstate) & exist([ls_spikes.sessionName '.behavior.mat']) & all(diff(intervals')>600)
        
    spkmat_ls = bz_SpktToSpkmat(ls_spikesBEHAV.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM_ls = bz_SpktToSpkmat(ls_spikesNREM.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM_post_ls = bz_SpktToSpkmat(ls_spikesNREM_post.times,'binSize',binSize(bins),'dt',binSize(bins));
    for spk = 1:size(spkmat_ls.data,2)
       spkmat_ls.zscoredData(:,spk) = zscore(spkmat_ls.data(:,spk));
       spkmatNREM_ls.zscoredData(:,spk) = zscore(spkmatNREM_ls.data(:,spk)); 
       spkmatNREM_post_ls.zscoredData(:,spk) = zscore(spkmatNREM_post_ls.data(:,spk)); 
    end
%     clear ls_spikes*

    

    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmat_ls.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmat_ls.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmat_ls.zscoredData(start:stop,:));
        if int == 1 
            [nah start] = min(abs(spkmatNREM_ls.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_ls.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_ls.zscoredData(start:stop,:),'rows','complete');
        elseif int == 3
            [nah start] = min(abs(spkmatNREM_post_ls.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_post_ls.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_post_ls.zscoredData(start:stop,:),'rows','complete');
        end
        corrmat{int}(abs(corrmat{int})>.99) = nan;
    end 
    clear spkmat*
    for i = 1:size(intervals,1)
        for j=1:size(intervals,1)
            temp = corrcoef(corrmat{i},corrmat{j},'rows','complete');
            Rmod(i,j) = temp(2);
        end
    end

    EV_ls(wind,bins,ii) = ((Rmod(2,3) - Rmod(2,1) * Rmod(1,3)) / sqrt((1 - Rmod(2,1)^2) * (1 - Rmod(1,3)^2)))^2;
    REV_ls(wind,bins,ii) = ((Rmod(2,1) - Rmod(2,3) * Rmod(1,3)) / sqrt((1 - Rmod(2,3)^2) * (1 - Rmod(1,3)^2)))^2;
    if EV_ls(wind,bins,ii) > 1
       error('problem..') 
    end
    end
    
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    hpcRegion{ii} = 'ca1';
    if isempty(hpc_spikes) % comment out to exclude CA3 recs
        hpc_spikes = bz_GetSpikes('region','ca3','noprompts',true);
        hpcRegion{ii} = 'ca3';
    end
    if ~isempty(hpc_spikes)  & ~isempty(it)
    hpc_spikesNREM = hpc_spikes;
    hpc_spikesBEHAV = hpc_spikes;

    for spk = 1:length(hpc_spikes.times)
%         if length(hpc_spikes.times{spk})./hpc_spikes.times{spk}(end) < 5 % 2.5Hz FR limit
       hpc_spikesNREM_post.times{spk} = Restrict(hpc_spikes.times{spk},it);   
%        hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},double(SleepState.ints.NREMstate));   
       hpc_spikesBEHAV.times{spk} = Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals); 
       hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
%         else
%        hpc_spikesNREM.times{spk} = [];   
%        hpc_spikesBEHAV.times{spk} = []; 
%        disp('excluded a cell...')
%         end
    end

    spkmat_hpc = bz_SpktToSpkmat(hpc_spikesBEHAV.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM_hpc = bz_SpktToSpkmat(hpc_spikesNREM.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM_post_hpc = bz_SpktToSpkmat(hpc_spikesNREM_post.times,'binSize',binSize(bins),'dt',binSize(bins));
    for spk = 1:size(spkmat_hpc.data,2)
    spkmat_hpc.zscoredData(:,spk) = zscore(spkmat_hpc.data(:,spk));
    spkmatNREM_hpc.zscoredData(:,spk) = zscore(spkmatNREM_hpc.data(:,spk)); 
    spkmatNREM_post_hpc.zscoredData(:,spk) = zscore(spkmatNREM_post_hpc.data(:,spk)); 
    end

    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmat_hpc.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmat_hpc.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmat_hpc.zscoredData(start:stop,:));
        if  int == 1 
            [nah start] = min(abs(spkmatNREM_hpc.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_hpc.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_hpc.zscoredData(start:stop,:),'rows','complete');
        elseif int == 3
            [nah start] = min(abs(spkmatNREM_post_hpc.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_post_hpc.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_post_hpc.zscoredData(start:stop,:),'rows','complete');            
        end
        corrmat{int}(abs(corrmat{int})>.99) = nan;
    end
    clear spkmat*
    for i = 1:size(intervals,1)
        for j=1:size(intervals,1)
            temp = corrcoef(corrmat{i},corrmat{j},'rows','complete');
            Rmod(i,j) = temp(2);
        end
    end

    EV_hpc(wind,bins,ii) = ((Rmod(2,3) - Rmod(2,1) * Rmod(1,3)) / sqrt((1 - Rmod(2,1)^2) * (1 - Rmod(1,3)^2)))^2;
    REV_hpc(wind,bins,ii) = ((Rmod(2,1) - Rmod(2,3) * Rmod(1,3)) / sqrt((1 - Rmod(2,3)^2) * (1 - Rmod(1,3)^2)))^2;
        
    if EV_hpc(wind,bins,ii) > 1
       error('problem..') 
    end
    end
    
    %% cross region
    spikes = bz_GetSpikes('noprompts',true);
    if ~isempty(spikes)  & ~isempty(it)
    spikesNREM = spikes;
    spikesBEHAV = spikes;
    for spk = 1:length(spikes.times)
%         spikesNREM.times{spk} = Restrict(spikes.times{spk},double(SleepState.ints.NREMstate));   
        spikesBEHAV.times{spk} = Restrict(spikes.times{spk},behavior.events.trialIntervals);  
        spikesNREM_post.times{spk} = Restrict(spikes.times{spk},it); 
        spikesNREM.times{spk} = Restrict(spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
       for spk2 = 1:length(spikes.times)
           if strcmp(spikes.region{spk},'hpc') & strcmp(spikes.region{spk2},'ls')
              idx(spk,spk2) = 1;
           else
               idx(spk,spk2) = 0;
           end
       end
    end

    spkmat = bz_SpktToSpkmat(spikesBEHAV.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM = bz_SpktToSpkmat(spikesNREM.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM_post = bz_SpktToSpkmat(spikesNREM_post.times,'binSize',binSize(bins),'dt',binSize(bins));
    for spk = 1:size(spkmat.data,2)
    spkmat.zscoredData(:,spk) = zscore(spkmat.data(:,spk));
    spkmatNREM.zscoredData(:,spk) = zscore(spkmatNREM.data(:,spk)); 
    spkmatNREM_post.zscoredData(:,spk) = zscore(spkmatNREM_post.data(:,spk)); 
    end
%     clear spikes*
    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmat.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmat.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmat.zscoredData(start:stop,:));
        if int == 1 
            [nah start] = min(abs(spkmatNREM.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM.zscoredData(start:stop,:),'rows','complete');
        elseif int == 3
            [nah start] = min(abs(spkmatNREM_post.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_post.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_post.zscoredData(start:stop,:),'rows','complete');
        end
        corrmat{int}(abs(corrmat{int})>.99) = nan;
        corrmat{int}(find(idx==0)) = nan;
    end
    clear idx

    for i = 1:size(intervals,1)
        for j=1:size(intervals,1)
            temp = corrcoef(corrmat{i},corrmat{j},'rows','complete');
            Rmod(i,j) = temp(2);
        end
    end

    EV_cross(wind,bins,ii) = ((Rmod(2,3) - Rmod(2,1) * Rmod(1,3)) / sqrt((1 - Rmod(2,1)^2) * (1 - Rmod(1,3)^2)))^2;
    REV_cross(wind,bins,ii) = ((Rmod(2,1) - Rmod(2,3) * Rmod(1,3)) / sqrt((1 - Rmod(2,3)^2) * (1 - Rmod(1,3)^2)))^2;

    end
    behavType{ii} = behavior.description;
    EV_ls(EV_ls==0) = nan;
    EV_hpc(EV_hpc==0) = nan;
    REV_ls(REV_ls==0) = nan;
    REV_hpc(REV_hpc==0) = nan;
    EV_cross(EV_cross==0) = nan;
    REV_cross(REV_cross==0) = nan;
    
    EV_ls(EV_ls==Inf) = nan;
    EV_hpc(EV_hpc==Inf) = nan;
    REV_ls(REV_ls==Inf) = nan;
    REV_hpc(REV_hpc==Inf) = nan;
    EV_cross(EV_cross==Inf) = nan;
    REV_cross(REV_cross==Inf) = nan;
    
    
%     figure(wind)
    
    subplot(3,4,bins)
    histogram(EV_ls(wind,bins,EV_ls(wind,bins,:)~=0)-REV_ls(wind,bins,EV_ls(wind,bins,:)~=0),-1:.1:1,'normalization','pdf','FaceColor','m')
    hold on
    histogram(EV_hpc(wind,bins,EV_hpc(wind,bins,:)~=0)-REV_hpc(wind,bins,EV_hpc(wind,bins,:)~=0),-1:.1:1,'normalization','pdf','FaceColor','k')
    hold off
    title(binSize(bins))
    subplot(3,4,10)
%     [cb] = cbrewer('qual','Set3',12,'pchip');
%     raincloud_plot('X',EV_cross(bins,:),'box_on', 1, 'bandwidth',.025,'color', [1 0 0], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .35+bins/2, 'dot_dodge_amount', .35+bins/2, 'box_col_match', 1,'cloud_edge_col', [1 0 0]);
% 
%     raincloud_plot('X',REV_cross(bins,:),'box_on', 1, 'bandwidth',.025,'color', [0 0 1], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .15+bins/2, 'dot_dodge_amount', .15+bins/2, 'box_col_match', 1,'cloud_edge_col', [0 0 1]);

%     boxplot(REV_ls(bins,:),'colors','k','symbol','')
%     hold on
%     plot(bins+.5,REV_cross(bins,:),'.b')
    boxplot(squeeze(EV_ls(wind,1:bins,:))','colors','b')
    boxplot(squeeze(EV_ls(wind,1:bins,:))','colors','b')
    
%     boxplot(binSize(1:bins),nanmedian(EV_cross(1:bins,:),2),sem(EV_cross(1:bins,:),2),'.b')
%     hold on
    errorbar(binSize(1:bins),nanmedian(REV_cross(wind,1:bins,:),3),sem((REV_cross(wind,1:bins,:)),3),'.r')
    title('hpc-ls cross region')
    set(gca,'xscale','log')
    hold off
    subplot(3,4,11)
%         raincloud_plot('X',EV_ls(bins,:),'box_on', 1, 'bandwidth',.025,'color', [1 0 0], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .35+bins/2, 'dot_dodge_amount', .35+bins/2, 'box_col_match', 1,'cloud_edge_col', [1 0 0]);
% 
%     raincloud_plot('X',REV_ls(bins,:),'box_on', 1, 'bandwidth',.025,'color', [0 0 1], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .15+bins/2, 'dot_dodge_amount', .15+bins/2, 'box_col_match', 1,'cloud_edge_col', [0 0 1]);

    errorbar(binSize(1:bins),nanmedian(EV_ls(wind,1:bins,:),3),sem((EV_ls(wind,1:bins,:)),3),'.b')
    hold on
    errorbar(binSize(1:bins),nanmedian(REV_ls(wind,1:bins,:),3),sem((REV_ls(wind,1:bins,:)),3),'.r')
    title('ls')
    set(gca,'xscale','log')
    hold off
    subplot(3,4,12)
%     raincloud_plot('X',EV_hpc(bins,:),'box_on', 1, 'bandwidth',.025,'color', [1 0 0], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .35+bins/2, 'dot_dodge_amount', .35+bins/2, 'box_col_match', 1,'cloud_edge_col', [1 0 0]);
% 
%     raincloud_plot('X',REV_hpc(bins,:),'box_on', 1, 'bandwidth',.025,'color', [0 0 1], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', .15+bins/2, 'dot_dodge_amount', .15+bins/2, 'box_col_match', 1,'cloud_edge_col', [0 0 1]);
    errorbar(binSize(1:bins),nanmedian(EV_hpc(wind,1:bins,:),3),sem((EV_hpc(wind,1:bins,:)),3),'.b')
    hold on
    errorbar(binSize(1:bins),nanmedian(REV_hpc(wind,1:bins,:),3),sem((REV_hpc(wind,1:bins,:)),3),'.r')
    hold off
    title('hpc')
    set(gca,'xscale','log')
    
    
    

    pause(.1)
    end
    wind = 1+wind;
    end
    end
cd ~/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
    end
