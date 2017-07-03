clear 
d = dir('*')
count = 1;

for i=3:length(d)
if d(i).isdir
    cd(d(i).name)
    xml = LoadParameters;
    lfp = bz_GetLFP([32:95]);
    rippleChan = bz_GetBestRippleChan(lfp);
    refChan = 3;
    lfp = bz_GetLFP([rippleChan refChan]);
    [b a] = butter(4,[130/(xml.lfpSampleRate./2) 180/(xml.lfpSampleRate./2)],'bandpass');
    for i=1:2
    lfp.filt(:,i) = FiltFiltM(b,a,double(lfp.data(:,i)));
    end

    % need to get intervals of SWS here...

    disp(['finding ripples with ch #: ' num2str(rippleChan) ', and noise ch #: ' num2str(refChan)])
    [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,'noise',...
        lfp.data(:,2),'durations',[12 100] ,'saveMat',false,'frequency',lfp.samplingRate);
    if length(ripples.times) > 2
    [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);

    ripples.rippleChan = rippleChan;
    ripples.refChan = refChan;

    save([xml.FileName '.hpc_ripples.event.mat'],'ripples','-v7.3')
    end
    
    cd ..
end
end