function [] = scratch_07_22_18_bayesianTemplate()
% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');
binSize = [.02];
tcFactors = [2:2:20];
count_ls = 1;
count_hpc = 1;
overlap = 5;



for ii=1
%     cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    spikes = bz_GetSpikes('noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    if exist([sessionInfo.FileName '.firingMaps.cellinfo.mat'])

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
    [firingMaps] = bz_firingMap1D(spikes,behavior,4);
%     binnedPhaseMaps = bz_phaseMap2Bins(phaseMaps.phaseMaps,firingMaps.rateMaps,behavior);
    
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
%        ls_spikesBEHAV.times{spk} = Restrict(ls_spikes.times{spk},behavior.events.trialIntervals);
       ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},[ripples.peaks-.2 ripples.peaks+.2]); 
    end
    end
   
    for bins = 1:length(binSize)
    if ~isempty(ls_spikes) & ~isempty(SleepState.ints.NREMstate) & exist([ls_spikes.sessionName '.behavior.mat']) & all(diff(intervals')>600)
        
%     spkmat_ls = bz_SpktToSpkmat(ls_spikesBEHAV.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
    spkmatNREM_ls = bz_SpktToSpkmat(ls_spikesNREM.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
    
    for t = 1:length(firingMaps.rateMaps)
%         template = squeeze(circ_mean(binnedPhaseMaps{t},[],2))+pi;
        template = squeeze(mean(firingMaps.rateMaps{t},2));
        for event = 1:length(ripples.peaks)
%             [n ts] = min(abs(spkmatNREM_ls.timestamps-ripples.peaks(event)));
            ts = round(ripples.peaks(event)*(1/spkmatNREM_ls.dt));
            if ts+20 < size(spkmatNREM_ls.data,1) & ts > 20
            data = spkmatNREM_ls.data(ts-20:ts+20,:);    
            [Pr_shuf prMax_shuf] = placeBayes((data')', bz_shuffleCircular(template(idx_ls,:)), spkmatNREM_ls.dt);
            [Pr, prMax] = placeBayes(data, template(idx_ls,:), spkmatNREM_ls.dt);
%             for i = 1:41
%             Pr(i,:) = minmax_norm(Pr(i,:));
%             Pr_shuf(i,:) = minmax_norm(Pr_shuf(i,:));
%             end
%             Pr_shuf= bz_shuffleCircular(Pr);  prMax_shuf = prMax(randperm(length(prMax)));
%             corrs_ls(t,event) = corr([1:41]',prMax,'rows','complete'); % exclude nans?
%             prMax(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0)=nan;
%             prMax_shuf(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0)=nan;
%             Pr(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0,:) = nan;
%             Pr_shuf(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0,:) = nan;
            if sum(sum(spkmatNREM_ls.data(ts-20:ts+20,:)))> 10 * overlap & sum(~isnan(sum(Pr')))>10
%                 corrs_ls_NaN(t,event) = corr([1:41]',prMax,'rows','complete');
%                 corrs_ls_NaN_shuf(t,event) = corr([1:41]',prMax_shuf,'rows','complete');
                [slope_ls(t,event) integral_ls(t,event)] = Pr2Radon(Pr);
                [slope_ls_shuf(t,event) integral_ls_shuf(t,event)] = Pr2Radon(Pr_shuf);
            else
%                 corrs_ls_NaN(t,event) = nan;
%                 corrs_ls_NaN_shuf(t,event) = nan;
                slope_ls(t,event) = nan;
                integral_ls(t,event) =nan;
                slope_ls_shuf(t,event) = nan;
                integral_ls_shuf(t,event) = nan;
            end
            
            else
%                 corrs_ls(t,event) = NaN;
%                 corrs_ls_NaN(t,event) = NaN;
%                 corrs_ls_NaN_shuf(t,event) = nan;
                slope_ls(t,event) = nan;
                integral_ls(t,event) =nan;
                slope_ls_shuf(t,event) = nan;
                integral_ls_shuf(t,event) = nan;
            
            end
        end
    end
%         preSleep_ls{count_ls} = corrs_ls(:,ripples.peaks<intervals(1,2));
%         postSleep_ls{count_ls} = corrs_ls(:,ripples.peaks>intervals(3,1));
%         behav_ls{count_ls} = corrs_ls(:,ripples.peaks>intervals(1,2) &...
%                              ripples.peaks<intervals(3,1));
%                          
%         preSleep_ls_NaN{count_ls} = corrs_ls_NaN(:,ripples.peaks<intervals(1,2));
%         postSleep_ls_NaN{count_ls} = corrs_ls_NaN(:,ripples.peaks>intervals(3,1));
%         behav_ls_NaN{count_ls} = corrs_ls_NaN(:,ripples.peaks>intervals(1,2) &...
%                              ripples.peaks<intervals(3,1));
%                          
%         preSleep_ls_NaN_shuf{count_ls} = corrs_ls_NaN_shuf(:,ripples.peaks<intervals(1,2));
%         postSleep_ls_NaN_shuf{count_ls} = corrs_ls_NaN_shuf(:,ripples.peaks>intervals(3,1));
%         behav_ls_NaN_shuf{count_ls} = corrs_ls_NaN_shuf(:,ripples.peaks>intervals(1,2) &...
%                              ripples.peaks<intervals(3,1));
                         
        ls_integral{count_ls} = integral_ls;
        ls_max_int{count_ls} = max(integral_ls);
        ls_max_int_shuf{count_ls} = max(integral_ls_shuf);
        ls_slope{count_ls} = slope_ls;
        ls_integral_shuf{count_ls} = integral_ls_shuf;
        ls_slope_shuf{count_ls} = slope_ls_shuf;
        count_ls = count_ls +1;       
        
        

%         subplot(3,2,3)
%         histogram(removeNAN(cell2vec(preSleep_ls_NaN)),'Normalization','pdf','FaceColor','b')
%         hold on
%         histogram(removeNAN(cell2vec(postSleep_ls_NaN)),'Normalization','pdf','FaceColor','r')
%         histogram(removeNAN(cell2vec(postSleep_ls_NaN_shuf)),'Normalization','pdf','FaceColor','k')
%         hold off
%         pause(.1) 
        clear corrs_ls* slope_ls* integral_ls*
    
    end
    
    
    
    if ~isempty(hpc_spikes) 
    lfp = bz_GetLFP(sessionInfo.ls);
    ls_power = zscore(fastrms(bz_Filter(double(lfp.data),'filter','butter','passband',[140 180],'order', 3),12));
    for event = 1:size(ripples.timestamps,1)
        start = round((ripples.peaks(event)-.025) * 1250);
        stop = round((ripples.peaks(event)+.025) * 1250);

        [ls_max{count_hpc}(event) a] = max(abs(ls_power(start:stop)));
    end
    
    hpc_spikesNREM = hpc_spikes;
%     hpc_spikesBEHAV = hpc_spikes;
    for spk = 1:length(hpc_spikes.times)
%         if length(hpc_spikes.times{spk})./hpc_spikes.times{spk}(end) < 5 % 2.5Hz FR limit
%        hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},double(SleepState.ints.NREMstate));   
%        hpc_spikesBEHAV.times{spk} = Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals); 
       hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},[ripples.peaks-.2 ripples.peaks+.2]); 
%         else
%        hpc_spikesNREM.times{spk} = [];   
%        hpc_spikesBEHAV.times{spk} = []; 
%        disp('excluded a cell...')
%         end
    end
%     spkmat_hpc = bz_SpktToSpkmat(hpc_spikesBEHAV.times,'overlap',overlap,'binSize',binSize(bins) * overlap);
    spkmatNREM_hpc = bz_SpktToSpkmat(hpc_spikesNREM.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
    for t = 1:length(firingMaps.rateMaps)
%         template = squeeze(circ_mean(binnedPhaseMaps{t},[],2))+pi;
        template = squeeze(mean(firingMaps.rateMaps{t},2));
        for event = 1:length(ripples.peaks)
%             [n ts] = min(abs(spkmatNREM_hpc.timestamps-ripples.peaks(event)));
            ts = round(ripples.peaks(event)*(1/spkmatNREM_hpc.dt));
            if ts+20 < size(spkmatNREM_hpc.data,1) & ts > 20
            data = spkmatNREM_hpc.data(ts-20:ts+20,:);    
            [Pr_shuf prMax_shuf] = placeBayes((data')', bz_shuffleCircular(template(idx_hpc,:)), spkmatNREM_hpc.dt);
            %shuffle template, not ripple, to preserve population burst?
            [Pr, prMax] = placeBayes(data, template(idx_hpc,:), spkmatNREM_hpc.dt);
%             for i = 1:41
%                 Pr(i,:) = minmax_norm(Pr(i,:));
%                 Pr_shuf(i,:) = minmax_norm(Pr_shuf(i,:));
%             end
%             Pr_shuf= bz_shuffleCircular(Pr); prMax_shuf = prMax(randperm(length(prMax)));
%             corrs_hpc(t,event) = corr([1:41]',prMax,'rows','complete'); % exclude nans?
%             prMax(sum(spkmatNREM_hpc.data(ts-20:ts+20,:)')==0)=nan;
%             prMax_shuf(sum(spkmatNREM_hpc.data(ts-20:ts+20,:)')==0)=nan;
%             Pr(sum(spkmatNREM_hpc.data(ts-20:ts+20,:)')==0,:) = nan;
%             Pr_shuf(sum(spkmatNREM_hpc.data(ts-20:ts+20,:)')==0,:) = nan;
            if sum(sum(spkmatNREM_hpc.data(ts-20:ts+20,:)))> 10 * overlap & sum(~isnan(sum(Pr')))>10
%                 corrs_hpc_NaN(t,event) = corr([1:41]',prMax,'rows','complete');
%                 corrs_hpc_NaN_shuf(t,event) = corr([1:41]',prMax_shuf,'rows','complete');
                  [slope_hpc(t,event) integral_hpc(t,event) ] = Pr2Radon(Pr);
                  [slope_hpc_shuf(t,event) integral_hpc_shuf(t,event)] = Pr2Radon(Pr_shuf);
            else 
%                 corrs_hpc_NaN(t,event) = nan;
%                 corrs_hpc_NaN_shuf(t,event) = nan;
                slope_hpc(t,event) = nan;
                integral_hpc(t,event) =nan;
                slope_hpc_shuf(t,event) = nan;
                integral_hpc_shuf(t,event) = nan;
            end
            else
%                 corrs_hpc(t,event) = NaN;
%                 corrs_hpc_NaN(t,event) = NaN;
%                 corrs_hpc_NaN_shuf(t,event) = nan;
                slope_hpc(t,event) = nan;
                integral_hpc(t,event) =nan;
                slope_hpc_shuf(t,event) = nan;
                integral_hpc_shuf(t,event) = nan;
            end
        end
    end
%         preSleep_hpc{count_hpc} = corrs_hpc(:,ripples.peaks<intervals(1,2));
%         postSleep_hpc{count_hpc} = corrs_hpc(:,ripples.peaks>intervals(3,1));
%         behav_hpc{count_hpc} = corrs_hpc(:,ripples.peaks>intervals(1,2) &...
%                              ripples.peaks<intervals(3,1));
%                          
%         preSleep_hpc_NaN{count_hpc} = corrs_hpc_NaN(:,ripples.peaks<intervals(1,2));
%         postSleep_hpc_NaN{count_hpc} = corrs_hpc_NaN(:,ripples.peaks>intervals(3,1));
%         behav_hpc_NaN{count_hpc} = corrs_hpc_NaN(:,ripples.peaks>intervals(1,2) &...
%                              ripples.peaks<intervals(3,1));
%                          
%         preSleep_hpc_NaN_shuf{count_hpc} = corrs_hpc_NaN_shuf(:,ripples.peaks<intervals(1,2));
%         postSleep_hpc_NaN_shuf{count_hpc} = corrs_hpc_NaN_shuf(:,ripples.peaks>intervals(3,1));
%         behav_hpc_NaN_shuf{count_hpc} = corrs_hpc_NaN_shuf(:,ripples.peaks>intervals(1,2) &...
%                              ripples.peaks<intervals(3,1));
                         
        hpc_integral{count_hpc} = integral_hpc;
        hpc_max_int{count_hpc} = max(integral_hpc);
        hpc_max_int_shuf{count_hpc} = max(integral_hpc_shuf);
        hpc_slope{count_hpc} = slope_hpc;
        hpc_integral_shuf{count_hpc} = integral_hpc_shuf;
        hpc_slope_shuf{count_hpc} = slope_hpc_shuf;
        count_hpc = count_hpc +1;       
        
        

%         subplot(3,2,4)
%         histogram(removeNAN(cell2vec(preSleep_hpc_NaN)),'Normalization','pdf','FaceColor','b')
%         hold on
%         histogram(removeNAN(cell2vec(postSleep_hpc_NaN)),'Normalization','pdf','FaceColor','r')
%         histogram(removeNAN(cell2vec(postSleep_hpc_NaN_shuf)),'Normalization','pdf','FaceColor','k')
%         hold off
%         pause(.1)
        clear corrs_hpc* slope_hpc* integral_hpc*
    end
    behavType{ii} = behavior.description;
   
    subplot(3,2,1)
    histogram(removeNAN(cell2vec(hpc_slope)),-.2:.01:.2,'Normalization','pdf','FaceColor','k')
    hold on
    histogram(removeNAN(cell2vec(ls_slope)),-.2:.01:.2,'Normalization','pdf','FaceColor','m')
    hold off
    title('slope')
    subplot(3,2,2)
    histogram(removeNAN(cell2vec(hpc_max_int)),'Normalization','pdf','FaceColor','k')
    hold on
    histogram(removeNAN(cell2vec(hpc_max_int_shuf)),'Normalization','pdf','FaceColor','r')
    hold off
    title('integral')
    subplot(3,2,5)
    histogram(removeNAN(cell2vec(ls_max_int)),'Normalization','pdf','FaceColor','b')
    hold on
    histogram(removeNAN(cell2vec(ls_max_int_shuf)),'Normalization','pdf','FaceColor','r')
    hold off
    title('integral')
    subplot(3,2,6)
    scatter(max(hpc_integral{count_hpc-1}),ls_max{count_hpc-1},'.k')
    hold on
    pause(.1)
    end
    end
    end
    
    save([sessionInfo.FileName '.bayesianResults.mat'],'-v7.3')
% cd ~/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
% savefig('/home/david/Dropbox/bayesianDecoder.fig')
% save('/home/david/Dropbox/bayesianDecoder_noExclusion.mat','-v7.3')
end
