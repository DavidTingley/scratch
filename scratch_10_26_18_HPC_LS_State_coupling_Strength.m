nrem_response = nan(98,10);
wake_response = nan(98,10);

d = dir('*201*');
for i=98:-1:1
cd(d(i).name)
popBursts = bz_LoadEvents(pwd,'popBursts');
ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
if ~isempty(ls_spikes) & ~isempty(popBursts) 
    if length(popBursts.bursts) > 100
SleepState = bz_LoadStates(pwd,'SleepState');
if isfield(SleepState.ints,'NREMstate')
spkMat = bz_SpktToSpkmat(ls_spikes,'binSize',.001,'overlap',1);
% save memory and add to single vector...
popRate = zeros(length(spkMat.data(:,1)),1);
for spk = 1:length(ls_spikes.times)
    popRate = popRate + zscore(Smooth(spkMat.data(:,spk),15));
end
popRate = popRate./spk;
for b = 1:length(popBursts.bursts)
ts = round(popBursts.bursts(b)*1000);
ls_resp(b) = max(popRate(ts-25:ts+25));
end
disc = (popBursts.amplitudes);
[nrem ] = InIntervals(popBursts.bursts,SleepState.ints.NREMstate);
[wake ] = InIntervals(popBursts.bursts,SleepState.ints.WAKEstate);
for dd= 1:10
    pct2 = prctile(popBursts.amplitudes,dd*10);
    pct1 = prctile(popBursts.amplitudes,(dd-1)*10);
idx = intersect(find(disc>pct1 & disc<pct2),find(nrem));
nrem_response(i,dd) = nanmean(ls_resp(idx));
idx = intersect(find(disc>pct1 & disc<pct2),find(wake));
wake_response(i,dd) = nanmean(ls_resp(idx));
end
% nrem_response(i,:) = zscore(nrem_response(i,:));
% wake_response(i,:) = zscore(wake_response(i,:));
clear ls_resp
subplot(3,2,1)
imagesc(nrem_response)
subplot(3,2,2)
imagesc(wake_response)
subplot(3,2,3)
cla
boundedline(1:size(nrem_response,2),nanmean(nrem_response),sem(nrem_response))
subplot(3,2,4)
cla
boundedline(1:size(nrem_response,2),nanmean(wake_response),sem(wake_response))
subplot(3,2,5)
plot(nanmean(nrem_response - wake_response))
pause(.1)
    end
    end
end
cd ~/datasets/ripples_LS/
end