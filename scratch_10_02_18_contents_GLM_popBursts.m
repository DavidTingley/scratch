function []  = scratch_10_02_18_contents_GLM_popBursts()
% d = dir('*201*');
count=0;
% [b a] = butter(4,[110/625 200/625],'bandpass');
[b a] = butter(4,[120/625 180/625],'bandpass');
    rec=1
%     
% for rec = 1:length(d)
%     cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    if  exist([sessionInfo.FileName '.placeFields.20_pctThresh.mat']) & exist([sessionInfo.FileName '.behavior.mat'])
        load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
%         load([sessionInfo.FileName '.olypherInfo_w_disc.cellinfo.mat'],'olypherInfo')
        load([sessionInfo.FileName '.rewardModulation.cellinfo.mat'],'rewardModulation')
%         load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
        SleepState = bz_LoadStates(pwd,'SleepState');
        load([sessionInfo.FileName '.behavior.mat'])
        for i=1:length(behavior.events.trials)
        bStart(i) = behavior.events.trials{i}.timestamps(1);
        bStop(i) = behavior.events.trials{i}.timestamps(end); 
        end
        pre = [0 min(bStart)];
        post = [max(bStop)  sessionInfo.recordingDuration];
        behav = [pre(end) post(1)]; clear bStop bStart
        
        % LOAD IN BAYESIAN DECODING SCORES HERE
        
        
        
        
        %% get initial correlations
        intervals = [pre; behav; post];
    
%         ca1 = load([sessionInfo.FileName '.CA1Ripples.events.mat']);
        popBursts = bz_LoadEvents(pwd,'popBursts');
        ca1.ripples = bz_LoadEvents(pwd,'CA1Ripples');
        if ~isempty(popBursts)
            ca1.ripples.peaks = popBursts.bursts;
            ca1.ripples.timestamps = popBursts.timestamps;
        end
        
        if ~isempty(SleepState) & ~isempty(SleepState.ints.NREMstate) & ~isempty(ca1.ripples) & isfield(SleepState.ints,'NREMstate')
            
                   
        spikes = bz_GetSpikes('noprompts',true);
        ls_spikes= bz_GetSpikes('noprompts',true,'region','ls');
        
        hpc_spikes= bz_GetSpikes('noprompts',true,'region','hpc');
        if isempty(hpc_spikes)
            hpc_spikes = bz_GetSpikes('noprompts',true,'region','ca3');
        end
        
        if ~isempty(ls_spikes) & ~isempty(hpc_spikes)
            

        % 
        hasField = nan(length(fields{1}),1);
        for i=1:length(fields)
            for j=1:length(fields{i})
                hasField(j) = hasField(j) + ~isempty(fields{i}{j});    
            end
        end
        %% get LFP power for all events..        
        ls.lfp = bz_GetLFP(sessionInfo.ls);
        if ~isempty(sessionInfo.ca1)
            ca1.lfp = bz_GetLFP(sessionInfo.ca1);
        else
            ca1.lfp = bz_GetLFP(sessionInfo.ca3);
        end
        ls_power = (fastrms(FiltFiltM(b,a,double(ls.lfp.data)),12));
        spkmat_hpc = bz_SpktToSpkmat(hpc_spikes,'binSize',.001,'overlap',1);
        spkmat_ls = bz_SpktToSpkmat(ls_spikes,'binSize',.001,'overlap',1);
        pr = zeros(size(spkmat_ls.data,1),1);
        for spk=1:size(spkmat_ls.data,2)
            pr = pr+zscore(Smooth(spkmat_ls.data(:,spk),10));
        end
        pr = pr./spk;
        
        pr_hpc = zeros(size(spkmat_hpc.data,1),1);
        for spk=1:size(spkmat_hpc.data,2)
            pr_hpc = pr_hpc+zscore(Smooth(spkmat_hpc.data(:,spk),10));
        end
        pr_hpc = pr_hpc./spk;
        
        
        clear s spkmat_* 
        hpc_power = (fastrms(FiltFiltM(b,a,double(ca1.lfp.data)),12));
        ls_rec = []; 
        hpc_rec = [];
%         hpc_pop = [];
        for event = 1:size(ca1.ripples.timestamps,1)
            start = round((ca1.ripples.timestamps(event,1)-.05) * 1250); % used to be 20 ms
            if start<1
                start = 1;
            end
            stop = round((ca1.ripples.timestamps(event,2)+.05) * 1250);
            if stop > length(ls_power)
                stop = length(ls_power);
            end
            [ls_max blah] = max(abs(ls_power(start:stop)));
            [hpc_max blah] = max(abs(hpc_power(start:stop)));
            hpc_rec=[hpc_rec;ls_max,hpc_max];
            
            
            start = round((ca1.ripples.timestamps(event,1)-.05) * 1000); % used to be 20 ms
            if start<1
                start = 1;
            end
            stop = round((ca1.ripples.timestamps(event,2)+.05) * 1000);
            if stop > length(ls_power)
                stop = length(ls_power);
            end
            [ls_pop(event)] = max(abs(pr(start:stop)));
            hpc_pop(event) =  max(abs(pr_hpc(start:stop)));
            
            if event > 1 & event < length(ca1.ripples.peaks)
                IRI(event) = (abs(ca1.ripples.peaks(event-1)-ca1.ripples.peaks(event)) + ...
                             abs(ca1.ripples.peaks(event+1)-ca1.ripples.peaks(event)))/2;
            else
                IRI(event) = nan;
            end
            
%             hold on
        end
        
        %% get that content, initialize
%         spatialContent = nan(length(spikes.times),length(ca1.ripples.peaks));
        rewardContent = nan(length(spikes.times),length(ca1.ripples.peaks));
        nSpikes = nan(length(spikes.times),length(ca1.ripples.peaks));
        PF = nan(length(spikes.times),length(ca1.ripples.peaks));
        cellLoc_wav = nan(length(spikes.times),length(ca1.ripples.peaks));
        
        % state stuff
        State = nan(length(ca1.ripples.peaks),1);
        SleepState.ints.NREMstate(:,2) = SleepState.ints.NREMstate(:,2) + 1;
        SleepState.ints.WAKEstate(:,2) = SleepState.ints.WAKEstate(:,2) + 1;
        SleepState.ints.REMstate(:,2) = SleepState.ints.REMstate(:,2) + 1;
        
        idx = InIntervals(ca1.ripples.peaks,SleepState.ints.NREMstate);
        State((idx)) = 1;
        idx = InIntervals(ca1.ripples.peaks,SleepState.ints.WAKEstate);
        State((idx)) = 2;
        idx = InIntervals(ca1.ripples.peaks,SleepState.ints.REMstate);
        State((idx)) = 3;
        
        % behavioral condition (pre/behav/sleep)
        condition = nan(length(ca1.ripples.peaks),1);
        for t = 1:3
        idx = InIntervals(ca1.ripples.peaks,intervals(t,:));
        condition(find(idx)) = t;
        end
        location = nan(length(ca1.ripples.peaks),1);
        for t = 1:3
        idx = InIntervals(ca1.ripples.peaks,intervals(t,:));
        if t == 2
            location(find(idx)) = 2;
        else
            location(find(idx)) = 1;
        end
        end
        
        
        % cell specific stuff
        for spk = 1:length(spikes.times)
%             if strcmp(hpc_spikes.region{spk},'hpc') | strcmp(hpc_spikes.region{spk},'ca3') | strcmp(hpc_spikes.region{spk},'ca1') 
%             rows = find(olypherInfo.results{spk}.discBins==2); 
%             cols = find(olypherInfo.results{spk}.smoothing==20);
%             cols = intersect(rows,cols);
            
            % olypher info variant
%             meanPeakRate(spk) = nanmax(olypherInfo.results{spk}.ratePeakInfo(cols)); % used to be nanmean
            
            for ind = 1:length(ca1.ripples.peaks)
                start = ((ca1.ripples.peaks(ind)-.05)); % used to be 20 ms
                if start<1
                    start = 1;
                end
                stop = ((ca1.ripples.peaks(ind)+.05));
                ripSpks = Restrict(spikes.times{spk},[start stop]);
                if ~isempty(ripSpks)
                    spkTimes(ind,spk) = mean(ripSpks);
                else
                    spkTimes(ind,spk) = nan;
                end
                duration(ind) = abs(diff(popBursts.timestamps(ind,:)));
%                 spatialContent(spk,ind) = meanPeakRate(spk);
                rewardContent(spk,ind) = rewardModulation.rewardGain(spk);
                PF(spk,ind) = (hasField(spk)>0);
                nSpikes(spk,ind) = length(ripSpks);
                if isfield(spikes,'chanDepthRelative_CA1PYR_wav')
                    cellLoc_wav(spk,ind) = spikes.chanDepthRelative_CA1PYR_wav(spk);
                end
                totSpks(spk,ind) = length(Restrict(spikes.times{spk},...
                                    [ca1.ripples.timestamps(ind,1)-.1 ca1.ripples.timestamps(ind,2)+.1]));
                
                
            end   
            
%             else 
%                 meanPeakRate = [];
%             end
        end
        
%         for ind = 1:length(ca1.ripples.peaks)
%            [a b] = sort(spkTimes(ind,:));
%            [a spikeOrder(ind,:)] = sort(b);
%         end
%         spikeOrder(isnan(spkTimes)) = 0;
%         u = unique(spikeOrder,'rows');
%         for uu=1:length(u)
%             for ind = 1:length(ca1.ripples.peaks)
%                 if sum(u(uu,:)==spikeOrder(ind,:))==size(spikeOrder,2)
%                     spkOrd(ind) = uu;
%                 end
% %                spkOrd(ind) = find(u==spikeOrder(ind,:)); 
%             end
%         end
        
    content.spikeTimes = spkTimes;
    content.interRipInterval{rec} = IRI;
    content.totalSpikes{rec} = totSpks;
    content.duration{rec} = duration;
    content.location{rec} = location;
    content.condition{rec} = condition; 
    content.SleepState{rec} = State;
    content.region{rec} = spikes.region;
    content.nCells{rec} = size(rewardContent,1);
%     content.spatialContent{rec} = spatialContent;
    content.rewardContent{rec} = rewardContent;
    content.nSpikes{rec} = nSpikes;
    content.PF{rec} = PF;
    content.hpc_popBurst{rec} = hpc_pop;
    content.ls_popBurst{rec} = ls_pop;
    content.cellLoc_wav{rec} = cellLoc_wav;
%     content.meanPeakRate{rec} = meanPeakRate;
    content.rewardGain{rec} = rewardModulation.rewardGain;
    content.hpc_power{rec} = hpc_rec(:,2);
    content.ls_power{rec} = hpc_rec(:,1);
    content.animal(rec) = sum(double(sessionInfo.animal));
    
    
%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% predictors = [nanmean(content.nSpikes{rec})' ...
%               nanmean(content.PF{rec})'...
%               nanmean(content.spatialContent{rec})' ... % check that this isn't just FR
%               nanmean(content.rewardContent{rec})'...
%               nanmean(content.cellLoc_wav{rec}.*content.nSpikes{rec})'...
%               content.condition{rec}...
%               content.location{rec}...
%               content.SleepState{rec}...
%               double(content.hpc_popBurst{rec})];
%           
% % actual = content.ls_popBurst{rec};
% actual = content.ls_power{rec}';
% 
% for p = 1:size(predictors,2)
%     if ~isnan(sum(predictors(:,p)))
%         [beta dev(p) ] = glmfit(predictors(:,p),actual,'normal');
%         yfit = glmval(beta,predictors(:,p),'identity');
%         sse(p) = sum((yfit-actual').^2);
% 
%         for iter = 1:100
%             pred_shuf = bz_shuffleCircular(predictors(:,p)')';
%     %         pred_shuf = predictors(randperm(size(predictors,1)),p);
%             [beta dev_shuf(p,iter)] = glmfit(pred_shuf,actual,'normal');
%             yfit = glmval(beta,pred_shuf,'identity');
%             sse_shuf(p,iter) = sum((yfit-actual').^2);
%         end
%     else
%         dev(p) = nan;
%         dev_shuf(p,1:100) = nan;
%     end
% end
% 
% 
% subplot(2,2,1)
% cla
% boundedline(1:size(predictors,2),nanmean(dev_shuf,2),std(dev_shuf,[],2))
% plot(dev,'r')
% subplot(2,2,2)
% devZ(rec,:) = (dev-nanmean(dev_shuf,2)')./std(dev_shuf,[],2)';
% imagesc(devZ)
% subplot(2,2,3)
% 
% subplot(2,2,4)
% cla
% boundedline(1:size(devZ,2),nanmean(devZ),nanstd(devZ))
% pause(.001)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([sessionInfo.FileName '.content_GLM_popBursts_50ms.mat'],'-v7.3')
    rec
    clear totSpks State meanPeakRate condition rewardModulation hpc_rec ls_rec meanPeakRate nSpikes *Content* PF cellLoc*
        end
        end
    end
%    cd /home/david/datasets/ripples_LS 
%    save('/home/david/datasets/ripples_LS/hpc_ripple_content_50ms_popBursts.mat','-v7.3')
% end

% for rec = 1:length(d)-1
%     if ~isempty(content.ls_power_z{rec})
% subplot(3,2,1)
% scatter(content.ls_power_z{rec},nanmean(content.PF{rec}),'.k')
% hold on
% ylabel('% of spks from place cells')
% xlabel('LS rippleband power')
% axis([0 40 0 .5])
% subplot(3,2,2)
% scatter(content.ls_power_z{rec},nanmean(content.rewardContent{rec}),'.k')
% hold on
% ylabel('reward gain probability')
% xlabel('LS rippleband power')
% % axis([0 40 0 .25])
% 
% 
% subplot(2,2,3)
% scatter(content.ls_power_z{rec},nanmean(content.nSpikes{rec}),'.k')
% hold on
%     end
% end