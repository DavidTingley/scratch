clear all
cd /home/david/datasets/ripples_LS
recordingList  = dir('*201*');
animals = {'DT1','DT2','DT5','DT7','DT8','DT9'};
[b a] = butter(4,[10/625 40/625],'bandpass');

winds = [10 25 40 60 80 100];

for window = 1:6
    
for i=1:length(recordingList)%hist:-1:45
   cd(recordingList(i).name) 
   sessionInfo = bz_getSessionInfo;
   ripples = bz_LoadEvents(pwd,'CA1Ripples');
   if ~isempty(ripples) & length(ripples.peaks) > 100
       ca1.channels = find(ismember(sessionInfo.region,{'hpc','ca1'}));
       ca1.lfp = bz_GetLFP(sessionInfo.channels(ca1.channels),'intervals',[ripples.peaks(:)-.10001 ripples.peaks(:)+.10001]);
       for rip = 1:length(ripples.peaks(:))
          for ch = 1:size(ca1.lfp(rip).data,2)
            [wavespec] = bz_WaveSpec(ca1.lfp(rip).data(:,ch),'samplingRate',1250,'frange',[120 200],'nfreqs',1,'ncyc',3,'space','lin');  
            ripPower(rip,ch,:) = abs(wavespec.data);
            
            [wavespec] = bz_WaveSpec(ca1.lfp(rip).data(:,ch),'samplingRate',1250,'frange',[10 40],'nfreqs',1,'ncyc',3,'space','lin');  
            sharpWavePower(rip,ch,:) = abs(wavespec.data);
            
            sharpWave(rip,ch,:) = FiltFiltM(b,a,double(ca1.lfp(rip).data(:,ch)));
          end
       end
       
       sessionInfo.chanDepthRelative_CA1PYR = nan(length(sessionInfo.channels),1);
   %%%
    for shank=1:length(sessionInfo.spikeGroups.groups)
%         subplot(4,3,shank)

        ind = find(sessionInfo.channels==sessionInfo.spikeGroups.groups{shank}(1));
        if ~isempty(ind)
        if strcmp(sessionInfo.region{ind},'hpc') | strcmp(sessionInfo.region{ind},'ca1')
        for chan=1:length(sessionInfo.spikeGroups.groups{shank})
            ind = find(ca1.lfp(rip).channels==sessionInfo.spikeGroups.groups{shank}(chan));
            if ~isempty(ind)
            [peakPower(chan) b] = max(squeeze(mean(ripPower(:,ind,:))));
            else
            peakPower(chan) = nan;
            end
        end
        if length(peakPower) > 4
%         hold on; plot(peakPower)
        % attempting ~1 um resolution...
%         try
        [blah peakChan] = max(smooth(interp1(peakPower,1:1/20:length(peakPower)),winds(window)));
%         plot(smooth(interp1(peakPower,1:1/20:length(peakPower)),55))
%         pause(.1)
%         [blah peakChan2] = max(smooth(interp1(peakPower(2:2:length(peakPower)),1:1/20:length(peakPower)/2),7));
%         catch
%         [blah peakChan] = max(smooth(interp(peakPower(1:2:length(peakPower))',5,1),5));
%         [blah peakChan2] = max(smooth(interp(peakPower(2:2:length(peakPower))',5,1),5));  
%         end
        
        peakChan = peakChan / 20;
%         peakChan2 = peakChan2  * 2 / 20;


%         peakChanLoc = find(ca1.lfp.channels==sessionInfo.spikeGroups.groups{shank}(peakChan));
        if peakChan >= 1.5 & peakChan <= length(peakPower) - 1.5 & blah > 800
           for chan=1:length(sessionInfo.spikeGroups.groups{shank})
                ind = find(sessionInfo.channels==sessionInfo.spikeGroups.groups{shank}(chan));
%                 if mod(chan,2)
                    dist = (peakChan - chan) * 10; % convert to um; + is above, - is below
                    sessionInfo.chanDepthRelative_CA1PYR(ind) = dist;
%                 else
%                     dist = (peakChan2 - chan) * 20; % convert to um; + is above, - is below
%                     sessionInfo.chanDepthRelative_CA1PYR(ind) = dist;
%                 end
            end
        else
            % we are likely not spanning CA1...
            for chan=1:length(sessionInfo.spikeGroups.groups{shank})
                ind = find(sessionInfo.channels==sessionInfo.spikeGroups.groups{shank}(chan));
                sessionInfo.chanDepthRelative_CA1PYR(ind) = nan;
            end
            
        end
        end
        clear peakPower
        end
        end
    end
    save([sessionInfo.FileName '.sessionInfo.mat'],'sessionInfo','-v7.3') 
       
   end
   sessionInfo.FileName
   cd /home/david/datasets/ripples_LS
end



d = dir('*201*');
for i=1:length(d)
cd(d(i).name)
sessionInfo = bz_getSessionInfo(pwd,'noprompts',true);
spikes = bz_GetSpikes('noprompt',true);
if ~isempty(spikes)
    if isfield(sessionInfo, 'chanDepthRelative_CA1PYR')
        spikes.chanDepthRelative_CA1PYR = nan(length(spikes.times),1);
        for spk = 1:length(spikes.times)
            if strcmp(spikes.region{spk},'ca1') | strcmp(spikes.region{spk},'hpc')
            chan = find(sessionInfo.channels ==spikes.maxWaveformCh(spk));
            spikes.chanDepthRelative_CA1PYR(spk) = sessionInfo.chanDepthRelative_CA1PYR(chan);
            end
        end
        save([sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes')
    end
end; cd /home/david/datasets/ripples_LS
end


                
depths =[];
depths_noWav =[];

for i=1:length(d)
cd(d(i).name)
spikes = bz_GetSpikes('noprompts',true); % requires buzcode
sessionInfo = bz_getSessionInfo;
if ~isempty(spikes) & isfield(spikes, 'chanDepthRelative_CA1PYR')
    for shank = 1:sessionInfo.spikeGroups.nGroups
        idx = sessionInfo.spikeGroups.groups{shank}+1;
        if ~isempty(dir([sessionInfo.FileName '.spk.' num2str(shank)]))
            if strcmp(sessionInfo.region{idx(1)},'ca1') | strcmp(sessionInfo.region{idx(1)},'hpc')
            try
            waves{shank} = LoadSpikeWaveforms([sessionInfo.FileName '.spk.' num2str(shank)],length(sessionInfo.spikeGroups.groups{shank}),32);
            clus{shank} = load([sessionInfo.FileName '.clu.' num2str(shank)]);
            catch 
            waves{shank} = [];
            clus{shank} = [];    
            end
        else
            waves{shank} = [];
            clus{shank} = [];
            end
        else
            waves{shank} = [];
            clus{shank} = [];
        end
    end
    for spk = 1:length(spikes.times)
        if strcmp(spikes.region{spk},'hpc') | strcmp(spikes.region{spk},'ca1')
            if ~isempty(waves{spikes.shankID(spk)})
            wav = waves{spikes.shankID(spk)};
            clu = clus{spikes.shankID(spk)}(2:end); 
            %% see if we have a channel above and below with the waveform on it...
            chan = spikes.maxWaveformCh(spk);
            ind = find(sessionInfo.spikeGroups.groups{spikes.shankID(spk)} == chan);
            spkOffset = 0;
            if ind ~= 1 & ind ~= length(sessionInfo.spikeGroups.groups{spikes.shankID(spk)})
                f = find(clu==spikes.cluID(spk));
                meanWaveforms = squeeze(mean(wav(f,:,:)));
                
                aboveChan = ind - 1;
%                 aboveChanID = find(sessionInfo.channels==sessionInfo.spikeGroups.groups{spikes.shankID(spk)}(aboveChan));
                belowChan = ind + 1;
%                 belowChanID = find(sessionInfo.channels==sessionInfo.spikeGroups.groups{spikes.shankID(spk)}(belowChan));

                aboveWaveform = (meanWaveforms(aboveChan,:));
                belowWaveform = (meanWaveforms(belowChan,:));
                

                [amp blah] = min(meanWaveforms(ind,:) - median(meanWaveforms(:)));
                [ampAbove blah] = min(aboveWaveform - median(meanWaveforms(:)));
                [ampBelow blah] = min(belowWaveform - median(meanWaveforms(:)));
                
                spkOffset = (ampAbove./amp - ampBelow./amp) * 10;
                
                if aboveChan ~= 1 & belowChan ~= size(meanWaveforms) % check one more for staggered probes
                    aboveChan = aboveChan - 1;
                    belowChan = belowChan + 1;
                    aboveWaveform = (meanWaveforms(aboveChan,:));
                    belowWaveform = (meanWaveforms(belowChan,:));
                    [ampAbove blah] = min(aboveWaveform - median(meanWaveforms(:)));
                    [ampBelow blah] = min(belowWaveform - median(meanWaveforms(:)));
                
                    spkOffset2 = (ampAbove./amp - ampBelow./amp) * 20;
                else
                    spkOffset2 = 0;
                end
                if abs(spkOffset) < abs(spkOffset2)
                    spkOffset = spkOffset2;
                end
            end
            %% add increment to channel depth
            depths = [depths;spikes.chanDepthRelative_CA1PYR(spk)+spkOffset];
            depths_noWav = [depths;spikes.chanDepthRelative_CA1PYR(spk)];
        end
    end
    end
    clear waves clus
    
end
cd /home/david/datasets/ripples_LS
end

figure(1)
subplot(3,3,window)
histogram(depths,-200:2:200)
title(winds(window))
figure(2)
subplot(3,3,window)
histogram(depths_noWav,-200:2:200)
title(winds(window))
pause(.1)
all_depths{window} = depths;
all_depths_noWav{window} = depths_noWav;
save('/home/david/Dropbox/depths_data.mat')
end
% 


