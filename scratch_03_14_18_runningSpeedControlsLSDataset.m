clear all
recordingList  = dir('*201*');


for i=1:length(recordingList)
   cd(recordingList(i).name) 
   sessionInfo = bz_getSessionInfo;
   if exist([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.behavior.mat'])
   load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
   [binnedPhaseMap] = bz_phaseMap2Bins(phaseMaps.phaseMaps,firingMaps.rateMaps,behavior);
   % smooth phases   
   for cond = 1:length(unique(behavior.events.trialConditions))
       for cell = 1:length(firingMaps.UID)
           for t = 1:length(find(behavior.events.trialConditions==cond))
               binnedPhaseMap{cond}(cell,t,isnan(binnedPhaseMap{cond}(cell,t,:)))=0;
               binnedPhaseMap{cond}(cell,t,:)=circ_smoothTS(squeeze(binnedPhaseMap{cond}(cell,t,:)),5,'method','mean','exclude',0);
           end
       end
   end
   % 
   for cond = 1:length(unique(behavior.events.trialConditions))
   if sum(behavior.events.trialConditions==cond) > 10
   trialIDX = find(behavior.events.trialConditions==cond);
   shuffleTrial = randperm(length(trialIDX));
   
   for cell = 1:length(firingMaps.UID)
       %% initialize
      corr_rate{i}{cond}(cell,1:length(trialIDX)) =nan;
      corr_delta_rate{i}{cond}(cell,1:length(trialIDX)) = nan;
      corr_phase{i}{cond}(cell,1:length(trialIDX)) = nan;
      corr_delta_phase{i}{cond}(cell,1:length(trialIDX)) = nan;
      % shuffled controls
      corr_rate_shuffle{i}{cond}(cell,1:length(trialIDX)) = nan;
      corr_delta_rate_shuffle{i}{cond}(cell,1:length(trialIDX)) = nan;
      corr_phase_shuffle{i}{cond}(cell,1:length(trialIDX)) = nan;
      corr_delta_phase_shuffle{i}{cond}(cell,1:length(trialIDX)) = nan;
      
   if strcmp(firingMaps.region{cell},'ls') && sum(sum(firingMaps.countMaps{cond}(cell,:,:))) > 1.5 * sum(behavior.events.trialConditions==cond)
   for t = 1:length(trialIDX)
      
      %% get velocity for all trials
      vel = abs(diff(behavior.events.trials{trialIDX(t)}.x)) + ...
          abs(diff(behavior.events.trials{trialIDX(t)}.y));
      
      %% get rate and change in rate
      rate = squeeze(firingMaps.rateMaps_box{cond}(cell,t,:));
      delta_rate = diff(rate);
      
      rate_shuffle = squeeze(firingMaps.rateMaps{cond}(cell,shuffleTrial(t),:));
      delta_rate_shuffle = diff(rate_shuffle);
      
      %% get phase and change in phase
      
      phase = squeeze(binnedPhaseMap{cond}(cell,t,:));
      delta_phase = diff(phase);
      
      phase_shuffle = squeeze(binnedPhaseMap{cond}(cell,shuffleTrial(t),:));
      delta_phase_shuffle = diff(phase_shuffle);
      
      %% calculate some things..
      vel = makeLength(vel,201);
      rate = makeLength(rate,201);
      delta_rate = makeLength(delta_rate,201);
      phase = makeLength(phase,201);
      delta_phase = makeLength(delta_phase,201);
      
      rate_shuffle = makeLength(rate_shuffle,201);
      delta_rate_shuffle = makeLength(delta_rate_shuffle,201);
      phase_shuffle = makeLength(phase_shuffle,201);
      delta_phase_shuffle = makeLength(delta_phase_shuffle,201);
      
      if sum(~isnan(phase)) > 5 & sum(phase~=0) > 5 & length(unique(phase)) > 2
      corr_rate{i}{cond}(cell,t) = corr(rate',vel');
      corr_delta_rate{i}{cond}(cell,t) = corr(delta_rate',vel');
      corr_phase{i}{cond}(cell,t) = circ_corrcl(phase',vel');
      if isinf(corr_phase{i}{cond}(cell,t))
          error
      end
      corr_delta_phase{i}{cond}(cell,t) = circ_corrcl(delta_phase',vel');
      % shuffled controls?
      
      corr_rate_shuffle{i}{cond}(cell,t) = corr(rate_shuffle',vel');
      corr_delta_rate_shuffle{i}{cond}(cell,t) = corr(delta_rate_shuffle',vel');
      corr_phase_shuffle{i}{cond}(cell,t) = circ_corrcl(phase_shuffle',vel');
      corr_delta_phase_shuffle{i}{cond}(cell,t) = circ_corrcl(delta_phase_shuffle',vel');
      end
   end
   end
   end
   end
   end
   end
   cd /home/david/datasets/lsDataset/
end


%% compile and plot
count = 1;
for i=1:length(recordingList)
    for cond = 1:length(corr_rate_shuffle{i})
        for cell = 1:size(corr_rate_shuffle{i}{cond},1)
            
        % means
        rate_corr(count) = nanmean(corr_rate{i}{cond}(cell,:));
        rate_corr_shuffle(count) = nanmean(corr_rate_shuffle{i}{cond}(cell,:));
        phase_corr(count) = nanmean(corr_phase{i}{cond}(cell,:));
        phase_corr_shuffle(count) = nanmean(corr_phase_shuffle{i}{cond}(cell,:));
        
        delta_rate_corr(count) = nanmean(corr_delta_rate{i}{cond}(cell,:));
        delta_rate_corr_shuffle(count) = nanmean(corr_delta_rate_shuffle{i}{cond}(cell,:));
        delta_phase_corr(count) = nanmean(corr_delta_phase{i}{cond}(cell,:));
        delta_phase_corr_shuffle(count) = nanmean(corr_delta_phase_shuffle{i}{cond}(cell,:));
        
        % stats
        [a pval_rate(count)] = ttest2(corr_rate{i}{cond}(cell,:),corr_rate_shuffle{i}{cond}(cell,:));
        [a pval_phase(count)] = ttest2(corr_phase{i}{cond}(cell,:),corr_phase_shuffle{i}{cond}(cell,:));
        [a pval_delta_rate(count)] = ttest2(corr_delta_rate{i}{cond}(cell,:),corr_delta_rate_shuffle{i}{cond}(cell,:));
        [a pval_delta_phase(count)] = ttest2(corr_delta_phase{i}{cond}(cell,:),corr_delta_phase_shuffle{i}{cond}(cell,:));
        count = 1+count;    
        end
    end
end
            
            
            
