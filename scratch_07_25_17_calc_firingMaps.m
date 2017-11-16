d  = dir('*201*');


 
%% compile data

for ii=1:length(d)
    cd(d(ii).name)
    if isempty(dir('*firingMaps*'))
        disp(['working on ' d(ii).name ''])
        
        load([d(ii).name '.behavior.mat'])
        load([d(ii).name '.sessionInfo.mat'])
        spikes = bz_GetSpikes('noprompts',true);
        if ~isempty(spikes)
        if ~isempty(sessionInfo.ca1)
            lfp = bz_GetLFP(sessionInfo.ca1);
        else
            lfp = bz_GetLFP(sessionInfo.ca3);
        end
        nBins = max(behavior.events.trials{1}.mapping);

        [firingMaps] = bz_firingMap1D(spikes,behavior,lfp,round(nBins./50),'savemat',true);
        for thresh = [.01 .05 .1 .2 .4]
            for cond = 1:length(firingMaps.rateMaps)
                fields{cond} = bz_getPlaceFields1D(firingMaps.rateMaps{cond},'minPeakRate',2);
            end
            save([d(ii).name '.placeFields.' num2str(thresh*100,'%0.2i') '_pctThresh.mat'],'fields')
        end
        end
%         end
   else
         spikes = bz_GetSpikes('noprompts',true);
        if ~isempty(spikes)
        load([d(ii).name '.firingMaps.cellinfo.mat'])
        for thresh = [.01 .05 .1 .2 .4]
            for cond = 1:length(firingMaps.rateMaps)
                fields{cond} = bz_getPlaceFields1D(firingMaps.rateMaps{cond},'minPeakRate',2,'percentThreshold',thresh);
            end
            save([d(ii).name '.placeFields.' num2str(thresh*100,'%0.2i') '_pctThresh.mat'],'fields')
        end
        end
    end
    cd /home/david/datasets/lsDataset
    clear fields firingMaps
end