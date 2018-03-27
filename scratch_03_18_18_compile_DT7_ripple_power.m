clear all
recCount = 1;
cd /home/david/datasets/ripples_LS
recordingList  = dir('*201*');

for i=1:length(recordingList)
   cd(recordingList(i).name) 
   sessionInfo = bz_getSessionInfo;
   if strcmp(sessionInfo.animal,'DT7')
       load([sessionInfo.FileName '.LSRipples.events.mat'])
%        states = bz_LoadStates(pwd,'SleepState');
%        nremTime = sum([100 diff(states.ints.NREMstate)']);
% bz_PlotRippleStats(ripples.maps,ripples.data,ripples.stats)
   if length(ripples.peaks) > 15
       avgWaveSpec{i} = zeros(length(ripples.peaks),64,100,625,'single');
       avgRippleBand{i} = zeros(length(ripples.peaks),64,625,'single');
       ripCount{i} = 0;
       
       for ripple = 1:length(ripples.peaks)
           
          lfp = bz_GetLFP(sessionInfo.spikeGroups.groups{7},'restrict',...
              [ripples.peaks(ripple)-.25 ripples.peaks(ripple)+.25]);
           for ch = 1:size(lfp.data,2)
               wavespec = bz_WaveSpec(lfp.data(:,ch),'frange',[1 200],'ncyc',3,...
                   'nfreqs',100,'space','lin','samplingRate',1250);
               avgWaveSpec{i}(ripple,ch,:,:) = abs(wavespec.data');
               avgRippleBand{i}(ripple,ch,:) = mean(abs(wavespec.data(:,60:100)'));
               ripCount{i} = ripCount{i} + 1;
%                imagesc(-.5:1/1250:.4999,sessionInfo.spikeGroups.groups{7},squeeze(mean(avgWaveSpec{i}(ripple,:,60:100,:),3)))
%                pause(.01)
           end
       end      
       numRips{i} = size(avgRippleBand{i},1);
       avgWaveSpec{i} =  squeeze(median(avgWaveSpec{i},1));
       avgRippleBand{i} =  squeeze(median(avgRippleBand{i},1));
       name = strsplit(sessionInfo.FileName,'_');
       name = strsplit(name{2},'u');
       depth(i) = str2num(name{1});
   bad = find(nanmean(avgRippleBand{i}')<70);
   for b = 1:length(bad)
      avgRippleBand{i}(bad(b),:) = nanmean(avgRippleBand{i}([bad(b)-1 bad(b)+1],:));   
   end
   subplot(5,4,recCount)
   imagesc(-.5:1/1250:.4999,sessionInfo.spikeGroups.groups{7},flipud(avgRippleBand{i}));
   
%    title(['day: ' num2str(recCount) ', depth: ' num2str(sessionInfo.depth) ', ripples: ' num2str(ripple)])
   recCount = 1+recCount
   end
   end
   pause(.1)
   cd /home/david/datasets/ripples_LS
end

save('/home/david/Dropbox/DT7_ripplePower_allChans.mat','-v7.3')