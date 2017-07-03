d = dir('*')
for i=17:length(d)
if d(i).isdir
    cd(d(i).name)
if ~exist([d(i).name '.hpc_RipplePhaseModulation.cellinfo.mat'])
% generateSessionInfo(s'DT7',0);
spikes = bz_GetSpikes('getwaveforms',true,'region','hpc','forcereload',true);
if ~isempty(spikes)
u = unique(spikes.maxWaveformCh);
for j=1:length(u)
lfp = bz_GetLFP(u(j));
f = find(spikes.maxWaveformCh==u(j));
t = bz_PhaseModulation(spikes.times(f),lfp,[120 180],'powerThresh',2,'plotting',false);
for k=1:length(f)
spikes.phaseModulation.phasedistros(:,f(k)) = t.phasedistros(:,k);
spikes.phaseModulation.phasebins = t.phasebins;
spikes.phaseModulation.spkphases{f(k)} = t.spkphases{k};
spikes.phaseModulation.detectorName = t.detectorName;
spikes.phaseModulation.detectorParams = t.detectorParams;
spikes.phaseModulation.phasestats.m(f(k))= t.phasestats.m(k);
spikes.phaseModulation.phasestats.r(f(k))= t.phasestats.r(k);
spikes.phaseModulation.phasestats.k(f(k))= t.phasestats.k(k);
spikes.phaseModulation.phasestats.p(f(k))= t.phasestats.p(k);
spikes.phaseModulation.phasestats.mode(f(k))= t.phasestats.mode(k);
end


save([spikes.sessionName '.hpc_RipplePhaseModulation.cellinfo.mat'],'spikes')
end
end
end
cd ..
clear spikes; i
end
end