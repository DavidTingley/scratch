d = dir('*201*');
   
for rec = 1:length(d)

   cd(d(rec).name)
   sessionInfo = bz_getSessionInfo;
   if exist([sessionInfo.FileName '.behavior.mat']) & exist([sessionInfo.FileName '.spikes.cellinfo.mat'])
       spikes = bz_GetSpikes;
       
       load([sessionInfo.FileName '.behavior.mat']);
       
       [times groups]=spikes2sorted(spikes.times);
       recDuration = (times(end));
       for spk = 1:length(spikes.times)
           for trial =1:length(behavior.events.trials)
              r(trial,spk) = length(Restrict(spikes.times{spk},[behavior.events.trials{trial}.timestamps(end) behavior.events.trials{trial}.timestamps(end)+1]));
              r_in(trial,spk) = length(Restrict(spikes.times{spk},[behavior.events.trials{trial}.timestamps(1) behavior.events.trials{trial}.timestamps(end)]));
              r_in(trial,spk) = r_in(trial,spk) ./ (behavior.events.trials{trial}.timestamps(end) - behavior.events.trials{trial}.timestamps(1));
           end
%            rewardGain(spk) = nanmean(r(:,spk))  ./  (length(spikes.times{spk})./recDuration);
           rewardGain(spk) = nanmean(r(:,spk))  ./  nanmean(r_in(:,spk));
       end
       rewardGain(isnan(rewardGain))=0;
       rewardGain(isinf(rewardGain))=1;

   rewardModulation.UID = spikes.UID;
   rewardModulation.region = spikes.region;
   rewardModulation.sessionname = spikes.sessionName;
   rewardModulation.rewardGain = rewardGain;
   save([spikes.sessionName '.rewardModulation.cellinfo.mat'])
   clear reward* r r_in recDuration times groups
   end
   cd /home/david/datasets/ripples_LS 
end