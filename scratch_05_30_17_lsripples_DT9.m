clear 
d = dir('*201*');
count = 1;

for i=1:length(d)
if d(i).isdir
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo;
    
    
    
    lfp = bz_GetLFP([sessionInfo.ca1 sessionInfo.refChan]);
%     ref = bz_GetLFP();
    
    
%     rippleChan = bz_GetBestRippleChan(lfp);
%     refChan = 3;
%     lfp = bz_GetLFP([rippleChan refChan]);

        [b a] = butter(4,[120/(xml.lfpSampleRate./2) 200/(xml.lfpSampleRate./2)],'bandpass');
    for ii=1:2
        lfp.filt(:,ii) = FiltFiltM(b,a,double(lfp.data(:,ii)));
    end

    % need to get intervals of SWS here...

%     disp(['finding ripples with ch #: ' num2str(rippleChan) ', and noise ch #: ' num2str(refChan)])
%     [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,'noise',...
%         lfp.data(:,2),'durations',[12 100],'saveMat',false,'frequency',lfp.samplingRate);
    
    [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,'noise',...
    lfp.data(:,2),'EMGThresh',.5,'thresholds',[2 4],'durations',[12 100],'saveMat',false,'frequency',lfp.samplingRate);

if length(ripples.times) > 2
    [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);

    ripples.rippleChan = rippleChan;
    ripples.refChan = refChan;

    save([sessionInfo.FileName '.CA1Ripples_4SD.events.mat'],'ripples','-v7.3')
    end
    
    cd ..
end
end