d = dir('*201*');
count=1;
count2=1;
% [b a] = butter(4,[110/625 200/625],'bandpass');
[b a] = butter(4,[120/625 180/625],'bandpass');
    
%     
for rec = length(d):-1:1
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
   spikes = bz_GetSpikes('noprompts',true);
   if ~isempty(spikes) & isfield(spikes,'chanDepthRelative_CA1PYR')
   for cell = 1:length(spikes.times)
       if ~isnan(spikes.chanDepthRelative_CA1PYR(cell)) &...
               strcmp(spikes.region{cell},'hpc') & ...
               isempty(sessionInfo.ca3) & length(spikes.times{cell})>500
       wavelet(count,:,:) = zeros(100,625);
%        if length(spikes.times{cell}) > 100
%            r = randperm(length(spikes.times{cell}));
%            r = r(1:100);
%        else
           r = randperm(length(spikes.times{cell}));
%        end
       
       lfp = bz_GetLFP(sessionInfo.ls,'intervals',...
           [spikes.times{cell}(r)-.25 spikes.times{cell}(r)+.25]); 
       for spk = 1:length(r)%length(spikes.times{cell})
           if length(lfp(spk).data) == 625
           wavespec = bz_WaveSpec(lfp(spk).data,'frange',[1 200],...
               'nfreqs',100,'ncyc',3,'space','lin','samplingRate',1250);
           wavelet(count,:,:) = squeeze(wavelet(count,:,:)) + abs(wavespec.data)';
           end
       end
       if isfield(spikes,'chanDepthRelative_CA1PYR_wav')
        depth(count) = spikes.chanDepthRelative_CA1PYR_wav(cell);
       else
        depth(count) = spikes.chanDepthRelative_CA1PYR(cell);   
       end
       for fr = 1:100
           wavelet_z(count,fr,:) = zscore(wavelet(count,fr,63:562));
       end

       
       count = 1+count;
       end
   end
   end
   
   
    %% CCG's here
    spikes = bz_GetSpikes('noprompt',true);
    if ~isempty(spikes)
    ls = bz_GetSpikes('noprompts',true,'region','ls');
    for spk = 1:length(spikes.times)
    if ~isempty(ls) & strcmp(spikes.region{spk},'hpc')
    [times groups] = spikes2sorted([{spikes.times{spk}},ls.times]);
    [ccg t] = CCG(times,groups,'binsize',.001,'duration',.5);
%     for k = 1:size(ccg,3)
    ccgs(count2,:) = squeeze(mean(ccg(:,1,:),3)); %used to sum over all pairs w/ pre-cell
    ccgs_norm(count2,:) = zscore(ccgs(count2,:));
    depths(count2) = spikes.chanDepthRelative_CA1PYR_wav(spk);
    count2 = 1+count2;
%     end
    end
    end
    end


    %% plotting here
    [a ord] = sort(depths);
    subplot(4,3,1)
    histogram(depths,50,'Normalization','pdf')
    subplot(4,3,2)
    plot(depths(ord),minmax_norm(smooth(mean(ccgs_norm(ord,255:265),2),200)))
    title('superficial cells impact LS more strongly')
    hold off
    if mod(count,10) == 0
    for s = 1:6
       subplot(4,3,s+2)
       p = prctile(depth,s*16.666);
       p1 = prctile(depth,(s-1)*16.666);
       idx = find(depth<=p & depth>p1);
       idx2 = find(depth>=p | depth<p1);
       
       imagesc(squeeze(mean(wavelet_z(idx,:,:),1))-squeeze(mean(wavelet_z(idx2,:,:),1)));
    %                for i = 1:100
    %                    temp(i,:) = zscore(temp(i,:));
    %                end
       caxis([-2 2])
       title(count)
       
       subplot(4,2,8)
       plot(squeeze(mean(mean(wavelet(idx,60:90,:),1))))
       hold on
        
    end
    hold off
    pause(.1)
    end
      
    

   cd E:\datasets\ripples_LS
    
end

timestamps = wavespec.timestamps;
freqs = wavespec.freqs;