


spikes = bz_GetSpikes;
isolationMetrics.UID = spikes.UID;
isolationMetrics.region = spikes.region;
isolationMetrics.sessionName = spikes.sessionName;

load([spikes.sessionName '.sessionInfo.mat'])

for i=1:length(sessionInfo.SpkGrps)
    f = find(spikes.shankID==i);
    [Fet, nFeatures] = LoadFeatures([sessionInfo.FileName '.fet.' num2str(i)]);
    nChannels = length(sessionInfo.SpkGrps(i).Channels);    
    nSamples = sessionInfo.SpkGrps(i).nSamples;    
    waveforms = LoadSpikeWaveforms([sessionInfo.FileName '.spk.' num2str(i)],nChannels,nSamples);
    clu = load([[sessionInfo.FileName '.clu.' num2str(i)]]);
    clu = clu(2:end);
    waveforms = single(reshape(waveforms,length(waveforms),nSamples*nChannels));
    
    for cell=1:length(f)
        isolationMetrics.isoDist_wav(f(cell)) = IsolationDistance(waveforms,find(clu==spikes.cluID(f(cell))));
        isolationMetrics.isoDist_fet(f(cell)) = IsolationDistance(Fet,find(clu==spikes.cluID(f(cell))));
        isolationMetrics.lratio_wav(f(cell)) = L_Ratio(waveforms,find(clu==spikes.cluID(f(cell))));
        isolationMetrics.lratio_fet(f(cell)) = L_Ratio(Fet,find(clu==spikes.cluID(f(cell))));
        isolationMetrics.waveformAmplitude(f(cell)) = max(abs(spikes.rawWaveform{f(cell)})) *.195;
    end   
    
end
isolationMetrics.waveformUnits = 'uV';
save([isolationMetrics.sessionName '.isolationMetrics.cellinfo.mat'],'isolationMetrics')
