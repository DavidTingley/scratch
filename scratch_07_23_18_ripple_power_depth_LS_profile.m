
d = dir('*201*');
for rec=83:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    spikes = bz_GetSpikes('noprompts',true);
    if ~isempty(spikes)
    spkmat = bz_SpktToSpkmat(spikes.times);
    ls_ripples = bz_LoadEvents(pwd,'LSRipples');
    ca1_ripples = bz_LoadEvents(pwd,'CA1Ripples');
    if ~isempty(ca1_ripples) & ~isempty(ls_ripples)
    for i=1:length(sessionInfo.region)
       if strcmp(sessionInfo.region{i},'ls')
           ch(i) = 1;
       else
           ch(i) = 0;
       end
    end
    idx = find(ch);
    lfp = bz_GetLFP(ls_ripples.rippleChan);
    power = zeros(length(ca1_ripples.timestamps),626);
    for rip = 1:length(ca1_ripples.timestamps)
        [a start] = min(abs(lfp.timestamps-(ca1_ripples.peaks(rip)-.25)));
        [a stop] = min(abs(lfp.timestamps-(ca1_ripples.peaks(rip)+.25)));
        while stop-start<625
           stop = stop+1; 
        end

       for ch = 1:size(lfp.data,2)
        wavespec = bz_WaveSpec(lfp.data(start:stop,ch),'samplingRate',1250,'nfreqs',1,'frange',[160 200],'space','lin');
        power(rip,:) = power(rip,:) + abs(wavespec.data)';
       end
       power(rip,:) = power(rip,:)./ch;
    end
    
      animal(rec) = sum(double(sessionInfo.animal));
      depth(rec) = sessionInfo.depth;
      power_all(rec,:) = mean(power);
      numRips(rec) = length(ca1_ripples.timestamps);
      recLength(rec) = max(spkmat.timestamps);
      noise(rec) = ca1_ripples.stdev;
      ca1_chan(rec) = ca1_ripples.rippleChan;
      ls_chan(rec) = ls_ripples.rippleChan;
      imagesc(power_all)
      pause(.1)
    end
    end
%     cd E:\datasets\ripples_LS
cd /home/david/datasets/ripples_LS
end

clf
u = unique(animal);
for i=2:length(u)
subplot(3,2,i-1)
f = find(animal == u(i));
plot(depth(f),power(f)./noise(f))
% plot(depth(f),numRips(f)./recLength(f))
end