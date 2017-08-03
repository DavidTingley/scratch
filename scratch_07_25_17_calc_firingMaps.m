d  = dir('*201*');


 
%% compile data

for ii=1:length(d)
    cd(d(ii).name)
    if isempty(dir('*firingMaps*'))
        disp(['working on ' d(ii).name ''])
        load([d(ii).name '.behavior.mat'])
        load([d(ii).name '.sessionInfo.mat'])
        spikes = bz_GetSpikes;
        lfp = bz_GetLFP(sessionInfo.thetaChans(2));
        nBins = max(behavior.events.trials{1}.mapping);

        [firingMaps] = bz_firingMap1D(spikes,behavior,lfp,round(nBins./50),'savemat',true);
    end
    cd /home/david/datasets/lsDataset
   
end