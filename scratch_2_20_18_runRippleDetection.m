d = dir('*201*');
count=1;
for rec = 1:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    sessionInfo.FileName 
    rippleChan = sessionInfo.ls;
    if ~isempty(rippleChan)
    refChan = sessionInfo.refChan;
    load([sessionInfo.FileName '.SleepState.states.mat'])
    
    lfp = bz_GetLFP([rippleChan refChan]);
    [b a] = butter(4,[120/(sessionInfo.lfpSampleRate./2) 180/(sessionInfo.lfpSampleRate./2)],'bandpass');
        for i=1:size(lfp.data,2)
            lfp.filt(:,i) = FiltFiltM(b,a,double(lfp.data(:,i)));
        end
%     if sum(diff(SleepState.ints.NREMstate')) > 60*3
        if ~isempty(refChan)
            [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,...
            'noise',lfp.data(:,2),...
            'durations',[10 100] ,...
            'passband', [120 180],...
            'threshold',[2 4],...
            'restrict',SleepState.ints.NREMstate,...
            'saveMat',false,...
            'frequency',lfp.samplingRate);
            if size(ripples.timestamps,1) > 1
                [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);
                ripples.rippleChan = rippleChan;
                ripples.refChan = refChan;

                save([sessionInfo.FileName '.LSRipples.events.mat'],'ripples','-v7.3')
            end
        else
            [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps, ...
            'durations',[10 100],...
            'passband', [120 180],...
            'threshold',[2 5],...
            'restrict',SleepState.ints.NREMstate,...
            'saveMat',false,...
            'frequency',lfp.samplingRate);
            if size(ripples.timestamps,1) > 1    
                [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);
                ripples.rippleChan = rippleChan;
                ripples.refChan = refChan;

                save([sessionInfo.FileName '.LSRipples.events.mat'],'ripples','-v7.3')
            end
        end
    else
        disp('not restricting, not enough sleep')
        if ~isempty(refChan)
            [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,...
            'noise',lfp.data(:,2),...
            'durations',[10 100] ,...
            'passband', [120 200],...
            'threshold',[2 4],...
            'saveMat',false,...
            'frequency',lfp.samplingRate);
            if size(ripples.timestamps,1) > 1
                [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);
                ripples.rippleChan = rippleChan;
                ripples.refChan = refChan;

                save([sessionInfo.FileName '.LSRipples.events.mat'],'ripples','-v7.3')
            end
        else
            [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps, ...
            'durations',[10 100],...
            'passband', [120 200],...
            'threshold',[2 4],...
            'saveMat',false,...
            'frequency',lfp.samplingRate);
            if size(ripples.timestamps,1) > 1    
                [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);
                ripples.rippleChan = rippleChan;
                ripples.refChan = refChan;

                save([sessionInfo.FileName '.LSRipples.events.mat'],'ripples','-v7.3')
            end
        end
    end
    end
   cd /home/david/datasets/ripples_LS 
end