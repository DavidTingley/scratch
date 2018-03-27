clear all
recCount = 1;
cd /home/david/datasets/ripples_LS
recordingList  = dir('*201*');
animals = {'DT1','DT2','DT5','DT7','DT8','DT9'};
colors = {'c','r','m','o','g','b'};

for i=1:length(recordingList)
   cd(recordingList(i).name) 
   sessionInfo = bz_getSessionInfo;
   if exist([sessionInfo.FileName '.LSRipples.events.mat'])
%    if strcmp(sessionInfo.animal,'DT7')
       load([sessionInfo.FileName '.LSRipples.events.mat'])
%        states = bz_LoadStates(pwd,'SleepState');
%        nremTime = sum([100 diff(states.ints.NREMstate)']);
   if length(ripples.peaks) > 15
       avgWaveSpec{i} = zeros(length(ripples.peaks),100,625,'single');
       avgRippleBand{i} = zeros(length(ripples.peaks),625,'single');
       ripCount{i} = 0;
       
       for ripple = 1:length(ripples.peaks)
          lfp = bz_GetLFP(sessionInfo.ls,'restrict',...
              [ripples.peaks(ripple)-.25001 ripples.peaks(ripple)+.25001]);
%            for ch = 1%1:size(lfp.data,2)
               wavespec = bz_WaveSpec(lfp.data,'frange',[1 200],'ncyc',3,...
                   'nfreqs',100,'space','lin','samplingRate',1250);
               avgWaveSpec{i}(ripple,:,:) = abs(wavespec.data');
               avgRippleBand{i}(ripple,:) = (mean(abs(wavespec.data(:,60:100)')));
               ripCount{i} = ripCount{i} + 1;
%            end
       end      
       avgWaveSpec{i} =  squeeze(median(avgWaveSpec{i},1));
       avgRippleBand{i} =  squeeze(median(avgRippleBand{i},1));
       avgRippleBand_z{i} =  zscore(squeeze(median(avgRippleBand{i},1)));
       
       depth(i) = sessionInfo.depth;
       animal{i} = sessionInfo.animal;    
       
       for ii = 1:6
           ind(ii) = strcmp(animal{i},animals{ii});
       end
       subplot(3,2,find(ind));
       
       plot(max(avgRippleBand{i}(300:325)),-depth(i),'.k')
       hold on
       pause(.001)
       
%    subplot(4,4,recCount)
%    imagesc(-.5:1/1250:.4999,sessionInfo.spikeGroups.groups{7},flipud(avgRippleBand{i}));
%    title(['day: ' num2str(recCount) ', depth: ' num2str(sessionInfo.depth) ', ripples: ' num2str(ripple)])
   recCount = 1+recCount
   end
   end
%    end
%    pause(.1)
   cd /home/david/datasets/ripples_LS
end

% save('/home/david/datasets/ripples_LS/DT7_ripplePower_allChans.mat','-v7.3')