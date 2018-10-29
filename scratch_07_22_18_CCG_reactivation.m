% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');

for int =1:3
    ccg_ls{int} = [];
    ccg_hpc{int} = [];
    ccg_cross{int} = [];
    hpc_rec{int} = [];
    placeFields{int} = [];
    ls_rec{int} = []; 
end

binSizes = [10 10 10];


for ii=1:length(d)
    cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    ls_spikesNREM_pre = ls_spikes;
    ls_spikesNREM_post = ls_spikes;
    ls_spikesBEHAV = ls_spikes;
    
    
    if exist([sessionInfo.FileName '.behavior.mat']) & ...
        ~isempty(ripples) & ...
        isfield(SleepState.ints,'NREMstate') & ~isempty(SleepState.ints.NREMstate) & ...
        exist([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
    load([sessionInfo.FileName '.behavior.mat'])
    load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
    pf = zeros(length(fields{1}),1);
    for c =1:length(fields)
        for sp = 1:length(fields{c})
            if ~isempty(fields{c}{sp})
                pf(sp) = pf(sp) + 1;
            end
        end
    end
    
    for i=1:length(behavior.events.trials)
       bStart(i) = behavior.events.trials{i}.timestamps(1);
       bStop(i) = behavior.events.trials{i}.timestamps(end); 
    end
    pre = [0 min(bStart)];
    post = [max(bStop)  sessionInfo.recordingDuration];
    behav = [pre(end) post(1)]; clear bStop bStart
    
    %% get initial correlations
    intervals = [pre; behav; post];
    
    % Skaggs 96 criteria for +/- 15 minutes
%     intervals(1,1) = intervals(1,2) - 60*15;
%     intervals(3,2) = intervals(3,1) + 60*15;

    
    if ~isempty(SleepState.ints.NREMstate) & ...
            exist([sessionInfo.FileName '.behavior.mat']) & all(diff(intervals')>300)
    if ~isempty(ls_spikes) 
        
        % ripples during NREM and intervals
%     for spk = 1:length(ls_spikes.times) 
%        ls_spikesBEHAV.times{spk} = [0;Restrict(ls_spikes.times{spk},behavior.events.trialIntervals)];
%        ls_spikesNREM_pre.times{spk} = [0;Restrict(Restrict(Restrict(ls_spikes.times{spk},[ripples.peaks-.2 ...
%            ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
%        ls_spikesNREM_post.times{spk} = [0;Restrict(Restrict(Restrict(ls_spikes.times{spk},[ripples.peaks-.2 ...
%            ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
%     end
    
%     only NREM sleep and intervals
     for spk = 1:length(ls_spikes.times)
       ls_spikesBEHAV.times{spk} = [spk;Restrict(ls_spikes.times{spk},behavior.events.trialIntervals)];
       ls_spikesNREM_pre.times{spk} = [spk;Restrict(Restrict((ls_spikes.times{spk} ...
           ),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
       ls_spikesNREM_post.times{spk} = [spk;Restrict(Restrict((ls_spikes.times{spk} ...
           ),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
     end

    % all states, just intervals
%     for spk = 1:length(ls_spikes.times)
%        ls_spikesBEHAV.times{spk} = [spk;Restrict(ls_spikes.times{spk},behavior.events.trialIntervals)];
%        ls_spikesNREM_pre.times{spk} = [spk;(Restrict(ls_spikes.times{spk} ...
%            ,intervals(1,:)))]; 
%        ls_spikesNREM_post.times{spk} = [spk;(Restrict(ls_spikes.times{spk} ...
%            ,intervals(3,:)))]; 
%     end


    for int = 1:size(intervals,1)
        if int == 1 
            [times groups]= spikes2sorted(ls_spikesNREM_pre.times);
        elseif int == 3
            [times groups]= spikes2sorted(ls_spikesNREM_post.times);
        elseif int == 2
            [times groups]= spikes2sorted(ls_spikesBEHAV.times);
        end
        [ccg{int} t] = CCG(times,groups,'binSize',.001,'duration',8);
    end
    for int = 1:size(intervals,1)
        for c1 = 1:size(ccg{int},2)
            for c2 = c1:size(ccg{int},3)
                if sum(ccg{1}(3801:4200,c1,c2)) > 2 & sum(ccg{2}(3801:4200,c1,c2)) > 2 & sum(ccg{3}(3801:4200,c1,c2)) > 2 & c1 ~= c2
                ccg_ls{int} = [ccg_ls{int}, ((squeeze(ccg{int}(:,c1,c2))))];
                ls_rec{int} = [ls_rec{int} ii];
                end
            end
        end
    end 
    clear ls_spikes* ccg
    end
        
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    hpcRegion{ii} = 'ca1';
    if isempty(hpc_spikes) % comment out to exclude CA3 recs
        hpc_spikes = bz_GetSpikes('region','ca3','noprompts',true);
        hpcRegion{ii} = 'ca3';
    end
    if ~isempty(hpc_spikes) 
%     for spk = 1:length(hpc_spikes.times)
%        hpc_spikesBEHAV.times{spk} = [0;Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals)];
%        hpc_spikesNREM_pre.times{spk} = [0;Restrict(Restrict(Restrict(hpc_spikes.times{spk},[ripples.peaks-.2 ...
%            ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
%        hpc_spikesNREM_post.times{spk} = [0;Restrict(Restrict(Restrict(hpc_spikes.times{spk},[ripples.peaks-.2 ...
%            ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
%     end
    
    % no ripple restriction
    for spk = 1:length(hpc_spikes.times)
       hpc_spikesBEHAV.times{spk} = [spk;Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals)];
       hpc_spikesNREM_pre.times{spk} = [spk;Restrict(Restrict((hpc_spikes.times{spk} ...
           ),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
       hpc_spikesNREM_post.times{spk} = [spk;Restrict(Restrict((hpc_spikes.times{spk} ...
           ),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
    end

      % all states, just intervals
%     for spk = 1:length(hpc_spikes.times)
%        hpc_spikesBEHAV.times{spk} = [spk;Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals)];
%        hpc_spikesNREM_pre.times{spk} = [spk;Restrict(hpc_spikes.times{spk} ...
%            ,intervals(1,:))]; 
%        hpc_spikesNREM_post.times{spk} = [spk;Restrict(hpc_spikes.times{spk} ...
%            ,intervals(3,:))]; 
%     end

    for int = 1:size(intervals,1)
        if int == 1 
            [times groups]= spikes2sorted(hpc_spikesNREM_pre.times);
        elseif int == 3
            [times groups]= spikes2sorted(hpc_spikesNREM_post.times);
        elseif int == 2
            [times groups]= spikes2sorted(hpc_spikesBEHAV.times);
        end
        [ccg{int} t] = CCG(times,groups,'binSize',.001,'duration',8);
    end
    for int = 1:size(intervals,1)
        for c1 = 1:size(ccg{int},2)
            for c2 = c1:size(ccg{int},3)
                if sum(ccg{1}(3801:4200,c1,c2)) > 2 & sum(ccg{2}(3801:4200,c1,c2)) > 2 & sum(ccg{3}(3801:4200,c1,c2)) > 2 & c1 ~= c2
                ccg_hpc{int} = [ccg_hpc{int}, ((squeeze(ccg{int}(:,c1,c2))))];
                placeFields{int} = [placeFields{int}; pf(c1) pf(c2)];
                hpc_rec{int} = [hpc_rec{int} ii];
                end
            end
        end
    end
    clear hpc_spikes* ccg
    end
    
    %% cross region
    spikes = bz_GetSpikes('noprompts',true);
    if ~isempty(spikes) 
        
%     for spk = 1:length(spikes.times)
%        spikesBEHAV.times{spk} = [0;Restrict(spikes.times{spk},behavior.events.trialIntervals)];
%        spikesNREM_pre.times{spk} = [0;Restrict(Restrict(Restrict(spikes.times{spk},[ripples.peaks-.2 ...
%            ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
%        spikesNREM_post.times{spk} = [0;Restrict(Restrict(Restrict(spikes.times{spk},[ripples.peaks-.2 ...
%            ripples.peaks+.2]),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
%     end

    % no ripple restriction
    for spk = 1:length(spikes.times)
       spikesBEHAV.times{spk} = [spk;Restrict(spikes.times{spk},behavior.events.trialIntervals)];
       spikesNREM_pre.times{spk} = [spk;Restrict(Restrict((spikes.times{spk} ...
           ),double(SleepState.ints.NREMstate)),intervals(1,:))]; 
       spikesNREM_post.times{spk} = [spk;Restrict(Restrict((spikes.times{spk} ...
           ),double(SleepState.ints.NREMstate)),intervals(3,:))]; 
    end

     % all states, just intervals
%       for spk = 1:length(spikes.times)
%        spikesBEHAV.times{spk} = [spk;Restrict(spikes.times{spk},behavior.events.trialIntervals)];
%        spikesNREM_pre.times{spk} = [spk;Restrict((spikes.times{spk} ...
%            ),intervals(1,:))]; 
%        spikesNREM_post.times{spk} = [spk;Restrict((spikes.times{spk} ...
%            ),intervals(3,:))]; 
%       end
    
    for int = 1:size(intervals,1)
        if int == 1 
            [times groups]= spikes2sorted(spikesNREM_pre.times);
        elseif int == 3
            [times groups]= spikes2sorted(spikesNREM_post.times);
        elseif int == 2
            [times groups]= spikes2sorted(spikesBEHAV.times);
        end
        [ccg{int} t] = CCG(times,groups,'binSize',.001,'duration',8);
    end
    for int = 1:size(intervals,1)
        for c1 = 1:size(ccg{int},2)
            for c2 = c1:size(ccg{int},3)
                if strcmp(spikes.region{c1},'ls') & sum(ccg{1}(3801:4200,c1,c2)) > 2 & sum(ccg{2}(3801:4200,c1,c2)) > 2 & sum(ccg{3}(3801:4200,c1,c2)) > 2 & c1 ~= c2
                    if strcmp(spikes.region{c2},'hpc') | strcmp(spikes.region{c2},'ca3')
                ccg_cross{int} = [ccg_cross{int}, ((squeeze(ccg{int}(:,c1,c2))))];
                    end
            end
        end
    end
    
    end
    clear spikes* ccg
    behavType{ii} = behavior.description;
      

    pause(.1)

    end
    end
    
    count = 1;
    warp = 2000; %[2000 3000 3500 3800 3900 3950];
    pre_hpc = nan(length(warp),1);
    post_hpc = nan(length(warp),1);
    pre_ls = nan(length(warp),1);
    post_ls = nan(length(warp),1);
    pre_hpc_shuf = nan(length(warp),1);
    post_hpc_shuf = nan(length(warp),1);
    pre_ls_shuf = nan(length(warp),1);
    post_ls_shuf = nan(length(warp),1);
    
% for w = 1:length(warp)
%     for c=1:size(ccg_hpc{1},2)
%         cc(:,c) = zscore(makeLength(Smooth(ccg_hpc{1}(warp(w):8001-warp(w),c),15),8001));
%         ccc(:,c) = zscore(makeLength(Smooth(ccg_hpc{3}(warp(w):8001-warp(w),c),15),8001));
%         cc_beh(:,c) = zscore(Smooth(ccg_hpc{2}(:,c),15));
%     end
%     pre_hpc(count) = corr2(cc_beh,cc);
%     post_hpc(count) = corr2(cc_beh,ccc);
%     pre_hpc_shuf(count) = corr2(bz_shuffleCircular(cc_beh),cc);
%     post_hpc_shuf(count) = corr2(bz_shuffleCircular(cc_beh),ccc);
%     subplot(7,2,1)
%     plot(4000./(4000-warp),pre_hpc,'.k')
%     hold on
%     plot(4000./(4000-warp),pre_hpc_shuf,'.r')
% %     hold off
% %     subplot(7,2,1)
%     plot(4000./(4000-warp),post_hpc,'+k')
%     hold on
%     plot(4000./(4000-warp),post_hpc_shuf,'+r')
%     xlabel('compression factor')
%     set(gca,'xscale','log')
%     hold off
%     title('HPC CCG reactivation')
%     subplot(7,2,3)
%     [a b o] = sort_cells(cc_beh',ccc',1);
%     imagesc(b)
%     subplot(7,2,5)
%     imagesc(a)
%     title('ccg HPC behav')
%     clear cc ccc cc_beh
%     
%     for c=1:size(ccg_ls{1},2)
%     cc(:,c) = zscore(makeLength(Smooth(ccg_ls{1}(warp(w):8001-warp(w),c),15),8001));
%     ccc(:,c) = zscore(makeLength(Smooth(ccg_ls{3}(warp(w):8001-warp(w),c),15),8001));
%     cc_beh(:,c) = zscore(Smooth(ccg_ls{2}(:,c),150));
%     end
%     pre_ls(count) = corr2(cc_beh,cc);
%     post_ls(count) = corr2(cc_beh,ccc);
%     pre_ls_shuf(count) = corr2(bz_shuffleCircular(cc_beh),cc);
%     post_ls_shuf(count) = corr2(bz_shuffleCircular(cc_beh),ccc);
%     subplot(7,2,2)
%     plot(4000./(4000-warp),pre_ls,'.m')
%     hold on
%     plot(4000./(4000-warp),pre_ls_shuf,'.r')
% %     hold off
% %     subplot(7,2,2)
%     plot(4000./(4000-warp),post_ls,'+m')
%     hold on
%     plot(4000./(4000-warp),post_ls_shuf,'+r')
%     set(gca,'xscale','log')
%     xlabel('compression factor')
%     hold off
%     title('LS CCG reactivation')
%     
%     subplot(7,2,4)
%     [a b o] = sort_cells(cc_beh',ccc',1);
%     imagesc(b)
%     subplot(7,2,6)
%     imagesc(a)
%     title('ccg LS behav')
%     
%     clear cc ccc cc_beh
%     
% 
%     count = 1 + count;
% end

% hpc_idx = hpc_rec{1} == ii;
hpc_pf_idx = sum(placeFields{1}'>0)>1;
hpc_idx = sum(placeFields{1}'>0)>1 & hpc_rec{1} == ii;
ls_idx = ls_rec{1} == ii;

if sum(hpc_idx(:))>2 & sum(ls_idx(:)) > 2
    subplot(7,2,7)
    scatter((mean(ccg_hpc{2}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_pf_idx))),...  % bias on track...
            (mean(ccg_hpc{1}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{1}(4001-200:4000,hpc_pf_idx))),'.k') % bias before
    xlabel('on track bias')
    ylabel('pre-behav bias')
    title(num2str(corr((mean(ccg_hpc{2}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_pf_idx)))',...  % bias on track...
            (mean(ccg_hpc{1}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{1}(4001-200:4000,hpc_pf_idx)))','rows','complete')))

    subplot(7,2,8)
    scatter((mean(ccg_hpc{2}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_pf_idx))),...  % bias on track...
            (mean(ccg_hpc{3}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{3}(4001-200:4000,hpc_pf_idx))),'.k') % bias before
    xlabel('on track bias')
    ylabel('post-behav bias')
    title(num2str(corr((mean(ccg_hpc{2}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_pf_idx)))',...  % bias on track...
            (mean(ccg_hpc{3}(4002:4001+200,hpc_pf_idx)) - mean(ccg_hpc{3}(4001-200:4000,hpc_pf_idx)))','rows','complete')))

    subplot(7,2,9)
    scatter((mean(ccg_ls{2}(4002:4001+200,:)) - mean(ccg_ls{2}(4001-200:4000,:))),...  % bias on track...
            (mean(ccg_ls{1}(4002:4001+200,:)) - mean(ccg_ls{1}(4001-200:4000,:))),'.m') % bias before
    xlabel('on track bias')
    ylabel('pre-behav bias')
    title(num2str(corr((mean(ccg_ls{2}(4002:4001+200,:)) - mean(ccg_ls{2}(4001-200:4000,:)))',...  % bias on track...
            (mean(ccg_ls{1}(4002:4001+200,:)) - mean(ccg_ls{1}(4001-200:4000,:)))','rows','complete')))

    subplot(7,2,10)
    scatter((mean(ccg_ls{2}(4002:4001+200,:)) - mean(ccg_ls{2}(4001-200:4000,:))),...  % bias on track...
            (mean(ccg_ls{3}(4002:4001+200,:)) - mean(ccg_ls{3}(4001-200:4000,:))),'.m') % bias before
    xlabel('on track bias')
    ylabel('post-behav bias')
    title(num2str(corr((mean(ccg_ls{2}(4002:4001+200,:)) - mean(ccg_ls{2}(4001-200:4000,:)))',...  % bias on track...
            (mean(ccg_ls{3}(4002:4001+200,:)) - mean(ccg_ls{3}(4001-200:4000,:)))','rows','complete')))

    subplot(7,2,11)
    N_p = union(find(mean(ccg_hpc{2}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_idx)) < 0 &...
                     mean(ccg_hpc{3}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{3}(4001-200:4000,hpc_idx)) < 0),...
                find(mean(ccg_hpc{2}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_idx)) > 0 &...
                     mean(ccg_hpc{3}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{3}(4001-200:4000,hpc_idx)) > 0));

    N_m = union(find(mean(ccg_hpc{2}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_idx)) < 0 &...
                     mean(ccg_hpc{3}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{3}(4001-200:4000,hpc_idx)) > 0),...
                find(mean(ccg_hpc{2}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{2}(4001-200:4000,hpc_idx)) > 0 &...
                     mean(ccg_hpc{3}(4002:4001+200,hpc_idx)) - mean(ccg_hpc{3}(4001-200:4000,hpc_idx)) < 0));

    plot(ii,length(N_p)./sum(hpc_idx),'.r')
    hold on
    plot(ii,length(N_m)./sum(hpc_idx),'.g')
    title('N++/N-- HPC')
    subplot(7,2,13)
    plot(ii,length(N_p),'.r')
    hold on
    plot(ii,length(N_m),'.g')
    title('N++/N-- HPC')
    
    subplot(7,2,12)
    N_p = union(find(mean(ccg_ls{2}(4001:4001+200,ls_idx)) - mean(ccg_ls{2}(4001-200:4001,ls_idx)) < 0 &...
                     mean(ccg_ls{3}(4001:4001+200,ls_idx)) - mean(ccg_ls{3}(4001-200:4001,ls_idx)) < 0),...
                find(mean(ccg_ls{2}(4001:4001+200,ls_idx)) - mean(ccg_ls{2}(4001-200:4001,ls_idx)) > 0 &...
                     mean(ccg_ls{3}(4001:4001+200,ls_idx)) - mean(ccg_ls{3}(4001-200:4001,ls_idx)) > 0));

    N_m = union(find(mean(ccg_ls{2}(4001:4001+200,ls_idx)) - mean(ccg_ls{2}(4001-200:4001,ls_idx)) < 0 &...
                     mean(ccg_ls{3}(4001:4001+200,ls_idx)) - mean(ccg_ls{3}(4001-200:4001,ls_idx)) > 0),...
                find(mean(ccg_ls{2}(4001:4001+200,ls_idx)) - mean(ccg_ls{2}(4001-200:4001,ls_idx)) > 0 &...
                     mean(ccg_ls{3}(4001:4001+200,ls_idx)) - mean(ccg_ls{3}(4001-200:4001,ls_idx)) < 0));
    plot(ii,length(N_p)./size(ccg_ls{1},2),'.r')
    hold on
    plot(ii,length(N_m)./size(ccg_ls{1},2),'.g')
    title('N++/N-- LS')
    
    subplot(7,2,14)
    plot(ii,length(N_p),'.r')
    hold on
    plot(ii,length(N_m),'.g')
    title('N++/N-- LS')


    pause(.01)
end
    end
cd ~/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
   
end
    
    
    
    
    
