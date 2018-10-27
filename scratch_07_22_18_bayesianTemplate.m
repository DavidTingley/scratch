function [] = scratch_07_22_18_bayesianTemplate()
% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');
binSize = [.01];
count_ls = 1;
count_hpc = 1;
overlap = 1;



% for ii=61:length(d)
%     cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    disp(['working on ' sessionInfo.FileName])
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    spikes = bz_GetSpikes('noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    popBursts = bz_LoadEvents(pwd,'popBursts');

    
    if exist([sessionInfo.FileName '.firingMaps.cellinfo.mat']) & ~isempty(popBursts)
%     ripples.peaks = popBursts.bursts;
%     ripples.timestamps = popBursts.timestamps;
    
    ls_spikesNREM = ls_spikes;
    ls_spikesBEHAV = ls_spikes;
    idx_ls = find(strcmp(spikes.region,'ls'));
    
    hpcRegion{count_hpc} = 'ca1';
    if isempty(hpc_spikes) % comment out to exclude CA3 recs
        hpc_spikes = bz_GetSpikes('region','ca3','noprompts',true);
        hpcRegion{count_hpc} = 'ca3';
    end
    idx_hpc = find(strcmp(spikes.region,'hpc') | strcmp(spikes.region,'ca3'));
    
    if exist([sessionInfo.FileName '.behavior.mat']) & ...
        ~isempty(ripples) & ...
        isfield(SleepState.ints,'NREMstate') & ~isempty(SleepState.ints.NREMstate)
    load([sessionInfo.FileName '.behavior.mat'])
    load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
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
    

    for bins = 1:length(binSize)

    
    
    if ~isempty(hpc_spikes) & ~isempty(ls_spikes)
    
    hpc_spikesNREM = hpc_spikes;
%     hpc_spikesBEHAV = hpc_spikes;
    for spk = 1:length(hpc_spikes.times)
       hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},[ripples.peaks-.5 ripples.peaks+.5]); 
    end
%     spkmat_hpc = bz_SpktToSpkmat(hpc_spikesBEHAV.times,'overlap',overlap,'binSize',binSize(bins) * overlap);
    spkmatNREM_hpc = bz_SpktToSpkmat(hpc_spikesNREM.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
    integral_hpc = nan(length(firingMaps.rateMaps),length(ripples.peaks));
    integral_hpc_shuf = nan(length(firingMaps.rateMaps),length(ripples.peaks),100);
    corrs_hpc = nan(length(firingMaps.rateMaps),length(ripples.peaks));
    corrs_hpc_shuf = nan(length(firingMaps.rateMaps),length(ripples.peaks),100);
    
    for t = 1:length(firingMaps.rateMaps)
        if size(firingMaps.rateMaps{t},2) >= 9
        rates = squeeze(mean(firingMaps.rateMaps{t}(idx_hpc,:,:),2));
        for spk = 1:size(rates,1)
            rates(spk,:) = mean_norm(rates(spk,:));
        end
        [a b ord] = sort_cells(rates);
%         
         for spk=1:length(spikes.times)
            pf(spk) = ~isempty(fields{t}{spk});
         end
         pf = find(pf);
         keep = intersect(pf,idx_hpc)-length(idx_ls); 
%          keep = 1:length(idx_hpc);
%          keep(sum(rates')==0) = [];
%         template = squeeze(circ_mean(binnedPhaseMaps{t},[],2))+pi;
%         template = squeeze(mean(firingMaps.rateMaps{t},2));
        
        template = squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_hpc,:,:),2));
        for i = 1:size(template,1)
           template(i,:) = mean_norm(Smooth(template(i,:),5)')';
        end
        for event = 1:length(ripples.peaks)
            % start bigger than the event itself
            start = round((round(ripples.timestamps(event,1) * 1000)-50) ./ (spkmatNREM_hpc.dt*1000));
            stop = round((round(ripples.timestamps(event,2) * 1000)+50) ./ (spkmatNREM_hpc.dt*1000));
            
            if stop < size(spkmatNREM_hpc.data,1) & stop-start > 5
                for spk = 1:size(spkmatNREM_hpc.data,2)
                        data(:,spk) = mean_norm(spkmatNREM_hpc.data(start:stop,spk)')';  
                        counts(:,spk) = (spkmatNREM_hpc.data(start:stop,spk)')';  
                end

                % get max rate bin
                [a b] = max(sum(data,2));
%                 or search from middle...
%                 b = round(size(data,1)/2);

                 % find clipping point at beginning
                 sta = 1;
                while sum(counts(b-sta,keep),2) > 1 & b-sta > 1 % max length is +/-50
                   sta = sta + 1;
                end
                 % find clipping point at beginning
                 sto = 1;
                while sum(counts(b+sto,keep),2) > 1 & sto +b < size(data,1)-1  % max length 500 ms
                   sto = sto + 1;
                end

                % redefine start/stop by pop burst..
                data = data(b-sta:b+sto,:);

                if size(data,1) > 5
                    [Pr, prMax] = placeBayes(data(:,keep), template(keep,:), spkmatNREM_hpc.dt*150);
    %                 % horse shit pfieffer/foster event clipping...
    %                 
    %                 gaps = unique([1 find(abs(diff(prMax))>15)' length(prMax)]);
    %                 [a b] = max(diff(gaps));
    %                 % only cuts if gap is found...
    %                 prMax = prMax(gaps(b):gaps(b+1));
    %                 data = data(gaps(b):gaps(b+1),:); 
                end
            if stop < size(spkmatNREM_hpc.data,1) & size(data,1) > 5
                % now calc correlations and control distro w/ cut data...
                corrs_hpc(t,event) = corr([1:length(prMax)]',prMax,'rows','complete');
                
                for iter = 1:100
%                     template_shuf = bz_shuffleCircular(squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_hpc,:,:),2)));
%                     for i = 1:size(template,1)
%                        template_shuf(i,:) = mean_norm(Smooth(template_shuf(i,:),5)')';
%                     end
                    template_shuf = bz_shuffleCellID(template);
                    [Pr_shuf prMax_shuf] = placeBayes((data(:,keep)')', template_shuf(keep,:), spkmatNREM_hpc.dt*150);
                    corrs_hpc_shuf(t,event,iter) = corr([1:length(prMax_shuf)]',prMax_shuf,'rows','complete');
                    subplot(4,2,6)
                    [slope_hpc_shuf(t,event,iter) integral_hpc_shuf(t,event,iter)] = Pr2Radon(Pr_shuf',1);
                end

                if sum(sum(spkmatNREM_hpc.data(start:stop,:)))> 5 * overlap & sum(~isnan(sum(Pr')))>5
                    subplot(4,2,4)
                    [slope_hpc(t,event) integral_hpc(t,event) ] = Pr2Radon(Pr',1);
                else 
                    corrs_hpc(t,event) = nan;
                    corrs_hpc_shuf(t,event,1:100) = nan;
                    slope_hpc(t,event) = nan;
                    integral_hpc(t,event) =nan;
                    slope_hpc_shuf(t,event) = nan;
                    integral_hpc_shuf(t,event,1:100) = nan;
                end
            else
                corrs_hpc(t,event) = nan;
                corrs_hpc_shuf(t,event,1:100) = nan;
                slope_hpc(t,event) = nan;
                integral_hpc(t,event) =nan;
                slope_hpc_shuf(t,event) = nan;
                integral_hpc_shuf(t,event,1:100) = nan;
                Pr = [];
            end
            
            nCells(event) = length(keep);
            spkCount(event) = sum(sum(data(:,keep)))./overlap;
            eventDuration(event) = (stop-start)*spkmatNREM_hpc.dt;
            end
            
            if ~isempty(Pr)
                
d = (integral_hpc-nanmean(integral_hpc_shuf,3))./nanstd(integral_hpc_shuf,[],3);

subplot(4,2,1)
histogram(integral_hpc,[0:.01:.3],'Normalization','pdf');
hold on
histogram(integral_hpc_shuf,[0:.01:.3],'Normalization','pdf')
title(['condition: ' num2str(t) ', event: ' num2str(event)])

% histogram(corrs_hpc,-1:.01:1,'Normalization','pdf');
% hold on
% histogram((corrs_hpc_shuf),-1:.01:1,'Normalization','pdf')
hold off

subplot(4,2,2)
% line([0 3],[.012 .012],'color','r');
bz_plotEphys([],'spikes',hpc_spikes,...
                    'scalelfp',5,...
                    'plotcells',hpc_spikes.UID(keep),...
                    'timewin',[ripples.timestamps(event,1)-.05 ripples.timestamps(event,2)+.05])

subplot(4,2,3)
scatter(max((d(:,1:event))),spkCount(1:event),'.k')
% scatter(absmax((corrs_hpc(:,1:event))),spkCount(1:event),'.k')
% xlabel('rank order corr')
ylabel('spk count')

subplot(4,2,4)
title([num2str(integral_hpc(t,event)) ' / ' num2str(d(t,event))])
% imagesc((Pr'))
% 
% title(corr(PR{count_hpc}(1:event)',max((d(:,1:event)))','rows','complete'))

subplot(4,2,5)
% plot(PR{count_hpc}(1:event),'.k')
plot(absmax((corrs_hpc(:,1:event))),'.k')
% imagesc(data(:,ord)')
title(corrs_hpc(t,event))

subplot(4,2,6)
title([num2str(nanmean(integral_hpc_shuf(t,event,:),3)) ' / ' num2str(d(t,event))])
% imagesc((Pr_shuf'))
% imagesc(rates(ord,:))
% nrem = InIntervals(ripples.peaks,SleepState.ints.NREMstate);
% wake = InIntervals(ripples.peaks,SleepState.ints.WAKEstate);
% errorbar(1,nanmean(max(integral_hpc(:,nrem))),nanstd(max(integral_hpc(:,nrem))'))
% hold on
% errorbar(2,nanmean(max(integral_hpc(:,wake))),nanstd(max(integral_hpc(:,wake))'))
% hold off        
% axis([0 3 0.004 .015])

subplot(4,2,7)
plot(max((d(:,1:event))),'.k')
% scatter(absmax((corrs_hpc(:,1:event))),popBursts.amplitudes(1:event),'.k')
xlabel('hpc replay')
ylabel('hpc rate')

subplot(4,2,8)
imagesc(flipud(data(:,keep)'))
% scatter(PR{count_hpc}(1:event),popBursts.amplitudes(1:event),'.k')
% xlabel('ls rate')
% ylabel('hpc rate')
pause(.001)
            end
clear data counts
        end
        end
    end
    
    
     behavType{count_hpc} = behavior.description;
   

    
    save([sessionInfo.FileName '.bayesianResults_ripples.mat'],'-v7.3')
    end
    end
    end
% cd ~/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
% savefig('/home/david/Dropbox/bayesianDecoder.fig')
% save('/home/david/Dropbox/bayesianDecoder_noExclusion.mat','-v7.3')

end
