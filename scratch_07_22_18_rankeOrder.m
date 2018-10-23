function [] = scratch_07_22_18_rankeOrder()
% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');
binSize = [.001];
count_ls = 1;
count_hpc = 1;
overlap = 1;

sessionInfo = bz_getSessionInfo;
% if ~exist(([sessionInfo.FileName '.rankOrder_ripples.mat']))

% for ii=61:length(d)
%     cd(d(ii).name)
  
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
    
    
%     if ~isempty(ls_spikes)
%     for spk = 1:length(ls_spikes.times)
% %        ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},double(SleepState.ints.NREMstate)); 
% %        ls_spikesBEHAV.times{spk} = Restrict(ls_spikes.times{spk},behavior.events.trialIntervals);
%        ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},[ripples.peaks-.25 ripples.peaks+.25]); 
%     end
%     end
%    for i=2:length(hpc_spikes.times)
    for bins = 1:length(binSize)
%     if ~isempty(ls_spikes) & ~isempty(SleepState.ints.NREMstate) & exist([ls_spikes.sessionName '.behavior.mat']) & all(diff(intervals')>600)
%         
% %     spkmat_ls = bz_SpktToSpkmat(ls_spikesBEHAV.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
%     spkmatNREM_ls = bz_SpktToSpkmat(ls_spikesNREM.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
%     
%     for t = 1:length(firingMaps.rateMaps)
% %         template = squeeze(circ_mean(binnedPhaseMaps{t},[],2))+pi;
% %         template = squeeze(mean(firingMaps.rateMaps{t},2));
%         template = squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_ls,:,:),2));
%         template_shuf = bz_shuffleCircular(squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_ls,:,:),2)));
%         for i = 1:size(template,1)
%            template(i,:) = Smooth(template(i,:),5);
%            template_shuf(i,:) = Smooth(template_shuf(i,:),5);
%         end
%         for event = 1:length(ripples.peaks)
% %             [n ts] = min(abs(spkmatNREM_ls.timestamps-ripples.peaks(event)));
%             ts = round(ripples.peaks(event)*(1/spkmatNREM_ls.dt));
%             if ts+20 < size(spkmatNREM_ls.data,1) & ts > 20
%             data = spkmatNREM_ls.data(ts-20:ts+20,:);    
% 
%             for iter = 1:100
%                 template_shuf = bz_shuffleCircular(squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_ls,:,:),2)));
%                 for i = 1:size(template,1)
%                    template_shuf(i,:) = Smooth(template_shuf(i,:),5);
%                 end
%                 [Pr_shuf prMax_shuf] = placeBayes((data')', template_shuf, spkmatNREM_ls.dt*5);
%                 [slope_ls_shuf(t,event,iter) integral_ls_shuf(t,event,iter)] = Pr2Radon(Pr_shuf);
%             end
%             
%             [Pr, prMax] = placeBayes(data, template, spkmatNREM_ls.dt*5);
% %             for i = 1:41
% %             Pr(i,:) = minmax_norm(Pr(i,:));
% %             Pr_shuf(i,:) = minmax_norm(Pr_shuf(i,:));
% %             end
% %             Pr_shuf= bz_shuffleCircular(Pr);  prMax_shuf = prMax(randperm(length(prMax)));
% %             corrs_ls(t,event) = corr([1:41]',prMax,'rows','complete'); % exclude nans?
% %             prMax(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0)=nan;
% %             prMax_shuf(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0)=nan;
% %             Pr(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0,:) = nan;
% %             Pr_shuf(sum(spkmatNREM_ls.data(ts-20:ts+20,:)')==0,:) = nan;
%             if sum(sum(spkmatNREM_ls.data(ts-20:ts+20,:)))> 10 * overlap & sum(~isnan(sum(Pr')))>10
% %                 corrs_ls_NaN(t,event) = corr([1:41]',prMax,'rows','complete');
% %                 corrs_ls_NaN_shuf(t,event) = corr([1:41]',prMax_shuf,'rows','complete');
%                 [slope_ls(t,event) integral_ls(t,event)] = Pr2Radon(Pr);
% %                 [slope_ls_shuf(t,event) integral_ls_shuf(t,event)] = Pr2Radon(Pr_shuf);
%             else
% %                 corrs_ls_NaN(t,event) = nan;
% %                 corrs_ls_NaN_shuf(t,event) = nan;
%                 slope_ls(t,event) = nan;
%                 integral_ls(t,event) =nan;
%                 slope_ls_shuf(t,event) = nan;
%                 integral_ls_shuf(t,event,1:100) = nan;
%             end
%             
%             else
% %                 corrs_ls(t,event) = NaN;
% %                 corrs_ls_NaN(t,event) = NaN;
% %                 corrs_ls_NaN_shuf(t,event) = nan;
%                 slope_ls(t,event) = nan;
%                 integral_ls(t,event) =nan;
%                 slope_ls_shuf(t,event) = nan;
%                 integral_ls_shuf(t,event,1:100) = nan;
%             
%             end
%         end
%     end
% %         preSleep_ls{count_ls} = corrs_ls(:,ripples.peaks<intervals(1,2));
% %         postSleep_ls{count_ls} = corrs_ls(:,ripples.peaks>intervals(3,1));
% %         behav_ls{count_ls} = corrs_ls(:,ripples.peaks>intervals(1,2) &...
% %                              ripples.peaks<intervals(3,1));
% %                          
% %         preSleep_ls_NaN{count_ls} = corrs_ls_NaN(:,ripples.peaks<intervals(1,2));
% %         postSleep_ls_NaN{count_ls} = corrs_ls_NaN(:,ripples.peaks>intervals(3,1));
% %         behav_ls_NaN{count_ls} = corrs_ls_NaN(:,ripples.peaks>intervals(1,2) &...
% %                              ripples.peaks<intervals(3,1));
% %                          
% %         preSleep_ls_NaN_shuf{count_ls} = corrs_ls_NaN_shuf(:,ripples.peaks<intervals(1,2));
% %         postSleep_ls_NaN_shuf{count_ls} = corrs_ls_NaN_shuf(:,ripples.peaks>intervals(3,1));
% %         behav_ls_NaN_shuf{count_ls} = corrs_ls_NaN_shuf(:,ripples.peaks>intervals(1,2) &...
% %                              ripples.peaks<intervals(3,1));
%                          
%         ls_rZ{count_ls} = (integral_ls - nanmean(integral_ls_shuf,3)) ./ nanstd(integral_ls_shuf,[],3);
%         ls_integral{count_ls} = integral_ls;
%         ls_max_int{count_ls} = max(integral_ls);
%         ls_max_int_shuf{count_ls} = max(integral_ls_shuf);
%         ls_slope{count_ls} = slope_ls;
%         ls_integral_shuf{count_ls} = integral_ls_shuf;
%         ls_slope_shuf{count_ls} = slope_ls_shuf;
%         count_ls = count_ls +1;       
%         
%         
% 
% %         subplot(4,2,3)
% %         histogram(removeNAN(cell2vec(preSleep_ls_NaN)),'Normalization','pdf','FaceColor','b')
% %         hold on
% %         histogram(removeNAN(cell2vec(postSleep_ls_NaN)),'Normalization','pdf','FaceColor','r')
% %         histogram(removeNAN(cell2vec(postSleep_ls_NaN_shuf)),'Normalization','pdf','FaceColor','k')
% %         hold off
% %         pause(.1) 
%         clear corrs_ls* slope_ls* integral_ls*
%     
%     end
    
    
    
    if ~isempty(hpc_spikes) & ~isempty(ls_spikes)
%     lfp = bz_GetLFP(sessionInfo.ls);
% %     hpc = bz_GetLFP(sessionInfo.ca1);
%     ls_power = zscore(fastrms(bz_Filter(double(lfp.data),'filter','butter','passband',[100 150],'order', 4),12));
%     
%     spkmat_ls = bz_SpktToSpkmat(ls_spikes,'binSize',.001,'overlap',1);
%     for i=1:size(spkmat_ls.data,2)
%         s(i,:) = zscore(Smooth(spkmat_ls.data(:,i),10));
%     end
%     pr = squeeze(mean(s));
% 
%     
%     for event = 1:size(ripples.peaks,1)
%         start = round((ripples.peaks(event)-.025) * 1250);
%         stop = round((ripples.peaks(event)+.025) * 1250);
%         
%         sta = round((ripples.peaks(event)-.025) * 1000);
%         sto = round((ripples.peaks(event)+.025) * 1000);
%         [PR{count_hpc}(event) a] = max(abs(pr(sta:sto)));
%         [ls_max{count_hpc}(event) a] = max(abs(ls_power(start:stop)));
%     end
    
    hpc_spikesNREM = hpc_spikes;
%     hpc_spikesBEHAV = hpc_spikes;
    for spk = 1:length(hpc_spikes.times)
%         if length(hpc_spikes.times{spk})./hpc_spikes.times{spk}(end) < 5 % 2.5Hz FR limit
%        hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},double(SleepState.ints.NREMstate));   
%        hpc_spikesBEHAV.times{spk} = Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals); 
       hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},[ripples.peaks-.5 ripples.peaks+.5]); 
%         else
%        hpc_spikesNREM.times{spk} = [];   
%        hpc_spikesBEHAV.times{spk} = []; 
%        disp('excluded a cell...')
%         end
    end
%     spkmat_hpc = bz_SpktToSpkmat(hpc_spikesBEHAV.times,'overlap',overlap,'binSize',binSize(bins) * overlap);
    spkmatNREM_hpc = bz_SpktToSpkmat(hpc_spikesNREM.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
    %% smooth ripple spike trains here
    for spk = 1:size(spkmatNREM_hpc.data,2)
        spkmatNREM_hpc.data(:,spk) = Smooth(spkmatNREM_hpc.data(:,spk),10);
    end
    
    rankOrder = nan(length(firingMaps.rateMaps),length(ripples.peaks));
    rankOrder_shuffle = nan(length(firingMaps.rateMaps),length(ripples.peaks),100);
    
   for t = 1:length(firingMaps.rateMaps)
        if size(firingMaps.rateMaps{t},2) >= 9
%         
         for spk=1:length(spikes.times)
            pf(spk) = ~isempty(fields{t}{spk});
         end
         pf = find(pf);
%          keep = 1:length(idx_hpc); 
         keep = intersect(pf,idx_hpc)-length(idx_ls); 
         
%         template = squeeze(circ_mean(binnedPhaseMaps{t},[],2))+pi;
%         template = squeeze(mean(firingMaps.rateMaps{t},2));
        
        template = squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_hpc,:,:),2));
        for i = 1:size(template,1)
           template(i,:) = mean_norm(Smooth(template(i,:),5)')';
        end
        
        for event = 1:length(ripples.peaks)
%             [n ts] = min(abs(spkmatNREM_hpc.timestamps-ripples.peaks(event)));
%             ts = round(ripples.peaks(event)*(1/spkmatNREM_hpc.dt));
%             back = ceil((ripples.peaks(event)-ripples.timestamps(event,1))./spkmatNREM_hpc.dt);
%             forward = abs(ceil((ripples.peaks(event)-ripples.timestamps(event,2))./spkmatNREM_hpc.dt));
            start = round(ripples.timestamps(event,1) * 1000)-50;
            stop = round(ripples.timestamps(event,2) * 1000)+50;
            
            if stop < size(spkmatNREM_hpc.data,1) & (stop-start) > 15
                for spk = 1:size(spkmatNREM_hpc.data,2)
                    data(:,spk) = mean_norm(spkmatNREM_hpc.data(start:stop,spk)')';  
                end 
                
                idx = intersect(find(nanmean(data)~=0),keep); % only take cells that spiked...
                idx2 = find(nanmean(data)~=0);
                
                [a b ord_template] = sort_cells(template(idx,:));
                [a b ord_data] = sort_cells(data(:,idx)');
                [a b ord_data_all] = sort_cells(data(:,idx)');

                if sum(sum(spkmatNREM_hpc.data(start:stop,:)))> 5 * overlap & length(idx) > 3                    
                    rankOrder(t,event) = corr(ord_template,ord_data,'rows','complete');
%                     rankOrder_allSPK(t,event) = corr(ord_template,ord_data_all,'rows','complete');
                    for iter = 1:100
                       template_shuf = bz_shuffleCircular(squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_hpc,:,:),2)));
                       for i = 1:size(template,1)
                          template_shuf(i,:) = mean_norm(Smooth(template_shuf(i,:),5)')';
                       end
%                        template_shuf = bz_shuffleCellID(template);
                       [a b ord_shuf] = sort_cells(template_shuf(idx,:));
                       rankOrder_shuffle(t,event,iter) = corr(ord_template,ord_shuf,'rows','complete');
                    end
                else 
                    rankOrder(t,event) = NaN;
                    rankOrder_shuffle(t,event,:) = nan(100,1);
                end
            else
                data = zeros(stop-start,size(spkmatNREM_hpc.data,2));
                rankOrder(t,event) = NaN;
                rankOrder_shuffle(t,event,:) = nan(100,1);
                idx = [];
            end
            nCells(t,event) = length(idx);
            spkCount(event) = sum(sum(spkmatNREM_hpc.data(start:stop,:)));
            eventDuration(event) = (stop-start)*spkmatNREM_hpc.dt;
%  
% subplot(4,2,1)
% histogram(rankOrder,[-1:.05:1],'Normalization','pdf');
% hold on
% histogram(rankOrder_shuffle,[-1:.05:1],'Normalization','pdf')
% hold off
% 
% subplot(4,2,2)
% % line([0 3],[.012 .012],'color','r');
% bz_MultiLFPPlot([],'spikes',hpc_spikes,...
%                     'spikeSpacingFactor',20,...  
%                     'scalelfp',5,...
%                     'timewin',[ripples.timestamps(event,:)],...
%                     'sortmetric',ord)
% 
% subplot(4,2,3)
% % plot(max((integral_hpc(:,1:event))),'.k');
% % imagesc(data')
% plot(absmax(rankOrder(:,1:event)),spkCount(1:event),'.k')
% xlabel('rank order corr')
% ylabel('spk count')
% 
% subplot(4,2,4)
% plot(absmax(rankOrder(:,1:event)),'.k')
% d = (rankOrder-mean(rankOrder_shuffle,3))./std(rankOrder_shuffle,[],3);
% % scatter(PR{count_hpc}(1:event),max((d(:,1:event))),'.k')
% % title(corr(PR{count_hpc}(1:event)',max((d(:,1:event)))','rows','complete'))
% 
% subplot(4,2,5)
% scatter(absmax(rankOrder(:,1:event)),eventDuration(1:event),'.k')
% xlabel('hpc rank order corr')
% ylabel('event duration (s)')
% 
% subplot(4,2,6)
% % scatter(rankOrder(t,1:event),rankOrder_allSPK(t,1:event),'.r')
% imagesc(rates(ord,:))
% % nrem = InIntervals(ripples.peaks,SleepState.ints.NREMstate);
% % wake = InIntervals(ripples.peaks,SleepState.ints.WAKEstate);
% % errorbar(1,nanmean(max(integral_hpc(:,nrem))),nanstd(max(integral_hpc(:,nrem))'))
% % hold on
% % errorbar(2,nanmean(max(integral_hpc(:,wake))),nanstd(max(integral_hpc(:,wake))'))
% % hold off        
% % axis([0 3 0.004 .015])
% 
% subplot(4,2,7)
% plot((rankOrder(1:t,1:event)),nCells(1:t,1:event),'.k')
% xlabel('rank order corr')
% ylabel('# of cells in event')
% 
% subplot(4,2,8)
% imagesc(data(:,ord(keep))')
% % scatter(PR{count_hpc}(1:event),popBursts.amplitudes(1:event),'.k')
% % xlabel('ls rate')
% % ylabel('hpc rate')
%  pause(.001)
clear data
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
%         hpc_rZ{count_hpc} = (integral_hpc - nanmean(integral_hpc_shuf,3)) ./ nanstd(integral_hpc_shuf,[],3);                 
%         hpc_integral{count_hpc} = integral_hpc;
%         hpc_max_int{count_hpc} = max(integral_hpc);
%         hpc_max_int_shuf{count_hpc} = max(integral_hpc_shuf);
%         hpc_slope{count_hpc} = slope_hpc;
%         hpc_integral_shuf{count_hpc} = integral_hpc_shuf;
%         hpc_slope_shuf{count_hpc} = slope_hpc_shuf;
        count_hpc = count_hpc +1;       
        
        

%         subplot(4,2,4)
%         histogram(removeNAN(cell2vec(preSleep_hpc_NaN)),'Normalization','pdf','FaceColor','b')
%         hold on
%         histogram(removeNAN(cell2vec(postSleep_hpc_NaN)),'Normalization','pdf','FaceColor','r')
%         histogram(removeNAN(cell2vec(postSleep_hpc_NaN_shuf)),'Normalization','pdf','FaceColor','k')
%         hold off
%         pause(.1)
%         clear corrs_hpc* slope_hpc* integral_hpc*
        
     behavType{count_hpc} = behavior.description;
   
%     subplot(4,2,1)
%     histogram(removeNAN(cell2vec(hpc_slope)),-.5:.01:.5,'Normalization','pdf','FaceColor','k')
%     hold on
%     histogram(removeNAN(cell2vec(ls_slope)),-.5:.01:.5,'Normalization','pdf','FaceColor','m')
%     hold off
%     title('slope')
%     subplot(4,2,2)
%     histogram(removeNAN(cell2vec(hpc_max_int)),'Normalization','pdf','FaceColor','k')
%     hold on
%     histogram(removeNAN(cell2vec(hpc_max_int_shuf)),'Normalization','pdf','FaceColor','r')
%     hold off
%     title('integral')
%     subplot(4,2,3)
%     histogram(removeNAN(cell2vec(hpc_max_int))-removeNAN(cell2vec(hpc_max_int_shuf)),'Normalization','pdf','FaceColor','k')
% %     subplot(4,2,4)
%     hold on
%     histogram(removeNAN(cell2vec(ls_max_int))-removeNAN(cell2vec(ls_max_int_shuf)),'Normalization','pdf','FaceColor','m')
%     hold off
%     subplot(4,2,5)
%     histogram(removeNAN(cell2vec(ls_max_int)),'Normalization','pdf','FaceColor','b')
%     hold on
%     histogram(removeNAN(cell2vec(ls_max_int_shuf)),'Normalization','pdf','FaceColor','r')
%     hold off
%     title('integral')
%     subplot(4,2,6)
%     scatter(max(hpc_rZ{count_hpc-1}),ls_max{count_hpc-1},'.k')
%     hold on
%     pause(.1)

    
    save([sessionInfo.FileName '.rankOrder_ripples.mat'],'-v7.3')
    end
    end
    end
    end
% cd ~/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
% savefig('/home/david/Dropbox/bayesianDecoder.fig')
% save('/home/david/Dropbox/bayesianDecoder_noExclusion.mat','-v7.3')
% end
end
