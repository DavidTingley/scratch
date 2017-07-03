d = dir('20170*')
count = 1;
for i=1:length(d)
whos
cd(d(i).name)
% generateSessionInfo('DT8',1);
% spikes = bz_GetSpikes('getwaveforms',false,'region','ls','forcereload',true);
if exist([d(i).name '.ls_RipplePhaseModulation.cellinfo.mat'])
    load([d(i).name '.ls_RipplePhaseModulation.cellinfo.mat'],'spikes')
    for j=1:length(spikes.times)
       num(count) = length(spikes.phaseModulation.spkphases{j});
       res(count) = spikes.phaseModulation.phasestats.r(j);
       m_angle(count) = spikes.phaseModulation.phasestats.m(j);
       p_val(count) = spikes.paseModulation.phasestats.p(j);
       shank(count) = spikes.shankID(j);
       clu(count) = spikes.cluID(j);
       wave(count,:) = spikes.rawWaveform{j};
       if num(count)>50
       polarplot(spikes.phaseModulation.phasestats.m(j),spikes.phaseModulation.phasestats.r(j),'.k')
       hold on
       end
       count = count + 1;
    end
end
cd ..
% clear spikes; i
end