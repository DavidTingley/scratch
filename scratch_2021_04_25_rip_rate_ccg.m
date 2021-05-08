d = dir('*201*');
for i=1:length(d)
   cd(d(i).name)
   ripples = bz_LoadEvents(pwd,'CA1Ripples');
   if exist([d(i).name '.spikes.cellinfo.mat']) & ~isempty(ripples)
   spikes = bz_GetSpikes('region','ca1','noPrompts',true);
   if isempty(spikes)
      spikes = bz_GetSpikes('region','hpc','noPrompts',true); 
   end
    if ~isempty(spikes)
%        [times groups] = spikes2sorted({ripples.peaks;spikes.spindices(:,1)});
%        [ccg t] = CCG(times,groups,'binSize',.01,'duration',60*60);
%        crossCo(i,:) = ccg(:,1,2);
        ripRate = hist(ripples.peaks,0:.1:ceil(spikes.spindices(end,1)));
        spkRate = hist(spikes.spindices(:,1),0:.1:ceil(spikes.spindices(end,1)));
        crossCo(i,:) = ccgBinned(ripRate,spkRate,10*60*60);
    end   
   end
   cd .. 
end