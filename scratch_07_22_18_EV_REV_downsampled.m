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
    for window = 1%0:300:3000
    % get the first XXX minutes of sleep..
%     temp = 1;
%     epoch = 0;
%     while sum(epoch) < 600 & temp < intervals(end)
%         it = SubtractIntervals(double(SleepState.ints.NREMstate),[0 intervals(2,2)+window; intervals(2,2)+temp+window intervals(end)]);
%         epoch = sum(diff(it'));
%         temp = temp+1;
%     end
    
    
    
    if ~isempty(ls_spikes) 
    for spk = 1:length(ls_spikes.times)
%        ls_spikesNREM_post.times{spk} = Restrict(ls_spikes.times{spk},it); 
       ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},double(SleepState.ints.NREMstate)); 
       ls_spikesBEHAV.times{spk} = Restrict(ls_spikes.times{spk},behavior.events.trialIntervals);
%        ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
    end
    end
   
    for bins = 1:length(binSize)
    if ~isempty(ls_spikes) & ~isempty(SleepState.ints.NREMstate) & exist([ls_spikes.sessionName '.behavior.mat']) & all(diff(intervals')>600)
        
        % downsampling to match spk counts
        for cell=1:length(ls_spikesNREM.times)
           preCount = length(Restrict(ls_spikesNREM.times{cell},intervals(1,:)));
           postCount = length(Restrict(ls_spikesNREM.times{cell},intervals(3,:))); 
           preRate = preCount ./ diff(intervals(1,:));
           postRate = postCount ./ diff(intervals(3,:));
           
           [a b] = max([preRate, postRate]);
           
           if b == 1  % pre > post
               if preCount > 0
               for s = 1:preCount
                   sub(s) = (preCount-s) ./ diff(intervals(1,:));
               end
               [a b] = min(abs(sub-postRate));
               r = randperm(preCount);
               preSpks = Restrict(ls_spikesNREM.times{cell},intervals(1,:));
               postSpks = Restrict(ls_spikesNREM.times{cell},intervals(3,:));
               preSpks(r(1:b))=[];
               end
           elseif b == 2 % pre < post
               if postCount > 0
               for s = 1:postCount
                   sub(s) = (postCount-s) ./ diff(intervals(3,:));
               end
               [a b] = min(abs(sub-preRate));
               r = randperm(postCount);
               preSpks = Restrict(ls_spikesNREM.times{cell},intervals(1,:));
               postSpks = Restrict(ls_spikesNREM.times{cell},intervals(3,:));
               postSpks(r(1:b))=[];
               end
           end
           ls_spikesNREM.times{cell} = [preSpks; postSpks];
           
           clear sub
        end        
        
    spkmatBEHAV_ls = bz_SpktToSpkmat(ls_spikesBEHAV.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM_ls = bz_SpktToSpkmat(ls_spikesNREM.times,'binSize',binSize(bins),'dt',binSize(bins));
%     spkmatNREM_post_ls = bz_SpktToSpkmat(ls_spikesNREM_post.times,'binSize',binSize(bins),'dt',binSize(bins));
    if ~isempty(spkmatNREM_ls.data)
    for spk = 1:size(spkmatBEHAV_ls.data,2)
       spkmatBEHAV_ls.zscoredData(:,spk) = zscore(spkmatBEHAV_ls.data(:,spk));
       spkmatNREM_ls.zscoredData(:,spk) = zscore(spkmatNREM_ls.data(:,spk)); 
%        spkmatNREM_post_ls.zscoredData(:,spk) = zscore(spkmatNREM_post_ls.data(:,spk)); 
    end
%     clear ls_spikes*

    

    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmatBEHAV_ls.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmatBEHAV_ls.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmatBEHAV_ls.zscoredData(start:stop,:));
        if int == 1 || int == 3
            [nah start] = min(abs(spkmatNREM_ls.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_ls.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_ls.zscoredData(start:stop,:),'rows','complete');
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
       warning('problem..') 
       EV_ls(wind,bins,ii) = nan;
       REV_ls(wind,bins,ii) = nan;
    end
    
    EV_ls(EV_ls==0) = nan;
    REV_ls(REV_ls==0) = nan;
    EV_ls(EV_ls==Inf) = nan;
    REV_ls(REV_ls==Inf) = nan;
    else
    EV_ls(wind,bins,ii) = nan; 
    REV_ls(wind,bins,ii) = nan;
    end
    end
    
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    hpcRegion{ii} = 'ca1';
    if isempty(hpc_spikes) % comment out to exclude CA3 recs
        hpc_spikes = bz_GetSpikes('region','ca3','noprompts',true);
        hpcRegion{ii} = 'ca3';
    end
    if ~isempty(hpc_spikes)  
   
        
    hpc_spikesNREM = hpc_spikes;
    hpc_spikesBEHAV = hpc_spikes;

    for spk = 1:length(hpc_spikes.times)
%         if length(hpc_spikes.times{spk})./hpc_spikes.times{spk}(end) < 5 % 2.5Hz FR limit

       hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},double(SleepState.ints.NREMstate));   
       hpc_spikesBEHAV.times{spk} = Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals); 
%        hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
%         else
%        hpc_spikesNREM.times{spk} = [];   
%        hpc_spikesBEHAV.times{spk} = []; 
%        disp('excluded a cell...')
%         end  
    end
    
    for cell=1:length(hpc_spikesNREM.times)
           preCount = length(Restrict(hpc_spikesNREM.times{cell},intervals(1,:)));
           postCount = length(Restrict(hpc_spikesNREM.times{cell},intervals(3,:))); 
           preRate = preCount ./ diff(intervals(1,:));
           postRate = postCount ./ diff(intervals(3,:));
           
           [a b] = max([preRate, postRate]);
           
           if b == 1  % pre > post
               if preCount > 0 
               for s = 1:preCount
                   sub(s) = (preCount-s) ./ diff(intervals(1,:));
               end
               [a b] = min(abs(sub-postRate));
               r = randperm(preCount);
               preSpks = Restrict(hpc_spikesNREM.times{cell},intervals(1,:));
               postSpks = Restrict(hpc_spikesNREM.times{cell},intervals(3,:));
               preSpks(r(1:b))=[];
               end
           elseif b == 2 % pre < post
               if postCount > 0
               for s = 1:postCount
                   sub(s) = (postCount-s) ./ diff(intervals(3,:));
               end
               [a b] = min(abs(sub-preRate));
               r = randperm(postCount);
               preSpks = Restrict(hpc_spikesNREM.times{cell},intervals(1,:));
               postSpks = Restrict(hpc_spikesNREM.times{cell},intervals(3,:));
               postSpks(r(1:b))=[];
               end
           end
           hpc_spikesNREM.times{cell} = [preSpks; postSpks];
           
           clear sub
    end 
        
    


    spkmatBEHAV_hpc = bz_SpktToSpkmat(hpc_spikesBEHAV.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM_hpc = bz_SpktToSpkmat(hpc_spikesNREM.times,'binSize',binSize(bins),'dt',binSize(bins));
    if ~isempty(spkmatNREM_hpc.data)
    for spk = 1:size(spkmatBEHAV_hpc.data,2)
        spkmatBEHAV_hpc.zscoredData(:,spk) = zscore(spkmatBEHAV_hpc.data(:,spk));
        spkmatNREM_hpc.zscoredData(:,spk) = zscore(spkmatNREM_hpc.data(:,spk)); 
    end

    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmatBEHAV_hpc.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmatBEHAV_hpc.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmatBEHAV_hpc.zscoredData(start:stop,:));
        if  int == 1 || int == 3
            [nah start] = min(abs(spkmatNREM_hpc.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_hpc.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_hpc.zscoredData(start:stop,:),'rows','complete');
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
       warning('problem..') 
       EV_hpc(wind,bins,ii) = nan;
       REV_hpc(wind,bins,ii) = nan;
       
    end
    EV_hpc(EV_hpc==0) = nan;
    REV_hpc(REV_hpc==0) = nan;
    EV_hpc(EV_hpc==Inf) = nan;
    REV_hpc(REV_hpc==Inf) = nan;
    
    else
    EV_hpc(wind,bins,ii) = nan; 
    REV_hpc(wind,bins,ii) = nan;
    end
    end
    
    %% cross region
    spikes = bz_GetSpikes('noprompts',true);
    if ~isempty(spikes) 
    spikesNREM = spikes;
    spikesBEHAV = spikes;
    
    for spk = 1:length(spikes.times)
    spikesNREM.times{spk} = Restrict(spikes.times{spk},double(SleepState.ints.NREMstate));   
    spikesBEHAV.times{spk} = Restrict(spikes.times{spk},behavior.events.trialIntervals);  
%         spikesNREM.times{spk} = Restrict(spikes.times{spk},[ripples.peaks-.1 ripples.peaks+.1]); 
   for spk2 = 1:length(spikes.times)
       if strcmp(spikes.region{spk},'hpc') & strcmp(spikes.region{spk2},'ls')
          idx(spk,spk2) = 1;
       else
           idx(spk,spk2) = 0;
       end
   end
    end
    
    for cell=1:length(spikesNREM.times)
           preCount = length(Restrict(spikesNREM.times{cell},intervals(1,:)));
           postCount = length(Restrict(spikesNREM.times{cell},intervals(3,:))); 
           preRate = preCount ./ diff(intervals(1,:));
           postRate = postCount ./ diff(intervals(3,:));
           
           [a b] = max([preRate, postRate]);
           
           if b == 1  % pre > post
               if preCount > 0
               for s = 1:preCount
                   sub(s) = (preCount-s) ./ diff(intervals(1,:));
               end
               [a b] = min(abs(sub-postRate));
               r = randperm(preCount);
               preSpks = Restrict(spikesNREM.times{cell},intervals(1,:));
               postSpks = Restrict(spikesNREM.times{cell},intervals(3,:));
               preSpks(r(1:b))=[];
               end
           elseif b == 2 % pre < post
               if postCount > 0
               for s = 1:postCount
                   sub(s) = (postCount-s) ./ diff(intervals(3,:));
               end
               [a b] = min(abs(sub-preRate));
               r = randperm(postCount);
               preSpks = Restrict(spikesNREM.times{cell},intervals(1,:));
               postSpks = Restrict(spikesNREM.times{cell},intervals(3,:));
               postSpks(r(1:b))=[];
               end
           end
           spikesNREM.times{cell} = [preSpks; postSpks];
           
           clear sub
    end 
        


    spkmatBEHAV = bz_SpktToSpkmat(spikesBEHAV.times,'binSize',binSize(bins),'dt',binSize(bins));
    spkmatNREM = bz_SpktToSpkmat(spikesNREM.times,'binSize',binSize(bins),'dt',binSize(bins));
    if ~isempty(spkmatNREM.data)
    for spk = 1:size(spkmatBEHAV.data,2)
    spkmatBEHAV.zscoredData(:,spk) = zscore(spkmatBEHAV.data(:,spk));
    spkmatNREM.zscoredData(:,spk) = zscore(spkmatNREM.data(:,spk)); 
    end
%     clear spikes*
    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmatBEHAV.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmatBEHAV.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmatBEHAV.zscoredData(start:stop,:));
        if int == 1 || int == 3
            [nah start] = min(abs(spkmatNREM.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM.zscoredData(start:stop,:),'rows','complete');
        end
        corrmat{int}(abs(corrmat{int})>.99) = nan;
        corrmat{int}(find(idx==0)) = nan;
    end
    

    for i = 1:size(intervals,1)
        for j=1:size(intervals,1)
            temp = corrcoef(corrmat{i},corrmat{j},'rows','complete');
            Rmod(i,j) = temp(2);
        end
    end

    EV_cross(wind,bins,ii) = ((Rmod(2,3) - Rmod(2,1) * Rmod(1,3)) / sqrt((1 - Rmod(2,1)^2) * (1 - Rmod(1,3)^2)))^2;
    REV_cross(wind,bins,ii) = ((Rmod(2,1) - Rmod(2,3) * Rmod(1,3)) / sqrt((1 - Rmod(2,3)^2) * (1 - Rmod(1,3)^2)))^2;
    EV_cross(EV_cross==0) = nan;
    REV_cross(REV_cross==0) = nan;
    EV_cross(EV_cross==Inf) = nan;
    REV_cross(REV_cross==Inf) = nan;
    
    else
    EV_cross(wind,bins,ii) = nan; 
    REV_cross(wind,bins,ii) = nan;
    end
    clear idx
    end
    
    behavType{ii} = behavior.description;
    
    
    
    
    
    
%     figure(wind)
    
    subplot(3,4,bins)
    histogram(EV_ls(wind,bins,EV_ls(wind,bins,:)~=0)-REV_ls(wind,bins,EV_ls(wind,bins,:)~=0),-1:.1:1,'normalization','pdf','FaceColor','m')
    hold on
    histogram(EV_hpc(wind,bins,EV_hpc(wind,bins,:)~=0)-REV_hpc(wind,bins,EV_hpc(wind,bins,:)~=0),-1:.1:1,'normalization','pdf','FaceColor','k')
    ylim([0 .5])
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
%     boxplot(squeeze(EV_ls(wind,1:bins,:))','colors','b')
%     boxplot(squeeze(EV_ls(wind,1:bins,:))','colors','b')
    
%     boxplot(binSize(1:bins),nanmedian(EV_cross(1:bins,:),2),sem(EV_cross(1:bins,:),2),'.b')
%     hold on
    errorbar(binSize(1:bins),nanmedian(EV_cross(wind,1:bins,:),3),sem((EV_cross(wind,1:bins,:)),3),'.b')
    hold on
    errorbar(binSize(1:bins),nanmedian(REV_cross(wind,1:bins,:),3),sem((REV_cross(wind,1:bins,:)),3),'.r')
    title('hpc-ls cross region')
    set(gca,'xscale','log')
    ylim([0 .5])
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
    ylim([0 .5])
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
    ylim([0 .5])
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
