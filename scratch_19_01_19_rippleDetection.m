cd /home/david/datasets/ripples_LS/
clear 
d = dir('*201*');
count = 1;

for i=1:length(d)
    if d(i).isdir
        cd(d(i).name)
        sessionInfo = bz_getSessionInfo;

        if ~isempty(sessionInfo.ca1) & ~isempty(sessionInfo.refChan)

        lfp = bz_GetLFP([sessionInfo.ca1 sessionInfo.refChan]);
    %     ref = bz_GetLFP();
        rips = bz_LoadEvents(pwd,'CA1Ripples');
if ~isempty(rips)
    %     rippleChan = bz_GetBestRippleChan(lfp);
    %     refChan = 3;
    %     lfp = bz_GetLFP([rippleChan refChan]);

            [b a] = butter(4,[120/(sessionInfo.lfpSampleRate./2) 200/(sessionInfo.lfpSampleRate./2)],'bandpass');
        for ii=1:2
            lfp.filt(:,ii) = FiltFiltM(b,a,double(lfp.data(:,ii)));
        end

        % need to get intervals of SWS here...

    %     disp(['finding ripples with ch #: ' num2str(rippleChan) ', and noise ch #: ' num2str(refChan)])
    %     [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,'noise',...
    %         lfp.data(:,2),'durations',[12 100],'saveMat',false,'frequency',lfp.samplingRate);
    if isfield(rips.detectorinfo.detectionparms,'EMGfilt')
        EMG = double(rips.detectorinfo.detectionparms.EMGfilt);
    elseif isfield(rips.detectorinfo.detectionparms,'EMGfilt')
        EMG = double(rips.detectorinfo.detectionparms.EMGThresh);
    else
        EMG = .75;
    end

        [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps, ...%'noise',lfp.data(:,2),...
                                   'EMGThresh',EMG,...
                                   'thresholds',[2 4],...
                                   'durations',[12 200],...
                                   'saveMat',false,...
                                   'restrict',rips.detectorinfo.detectionparms.restrict,...
                                   'frequency',lfp.samplingRate);

        if length(ripples.timestamps) > 2
            [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);

            ripples.rippleChan = sessionInfo.ca1;
            ripples.refChan = sessionInfo.refChan;

            save([sessionInfo.FileName '.CA1Ripples_4SD.events.mat'],'ripples','-v7.3')
        end
end
        end
        cd /home/david/datasets/ripples_LS/
    end
end