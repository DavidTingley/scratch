d = dir('*201*');
b = dir('*/*behavior*');
l = 1; h=1;
tau = 5;

for i=1:length(d)
%     load([d(i).folder '/' d(i).name]);
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo;

    if exist([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
    b = dir('*behavior*');
    load([b(1).name]);
    load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
    load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
    load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
    spikes = bz_GetSpikes;
%     cd /home/david/datasets/lsDataset
%     if strcmp(behavior.description,'wheel alternation')
    [firingMaps] = bz_firingMap1D(spikes,behavior,tau);
        [binnedfiringMaps.phaseMaps] = bz_phaseMap2Bins(phaseMaps.phaseMaps,firingMaps.rateMaps,behavior);
    for cell = 1:length(firingMaps.UID)
        for condition = 1:length(firingMaps.rateMaps)
            if sum(sum(firingMaps.countMaps{condition}(cell,:,:))) > 1.5 * size(firingMaps.countMaps{condition},2)
           if strcmp(firingMaps.region{cell},'ls') & size(firingMaps.rateMaps{condition},2) > 10
%             l_maps(l,:) = (squeeze(mean(firingMaps.rateMaps{condition}(cell,:,:),2)));
            ll = (squeeze((squeeze(binnedfiringMaps.phaseMaps{condition}(cell,:,:)))));
            for trial = 1:size(ll,1)
               ll(trial,:) = circ_smoothTS(ll(trial,:),tau,'method','mean','exclude',0); 
            end
            l_rate_maps_smooth(l,:) = (squeeze(mean(firingMaps.rateMaps_box{condition}(cell,:,:),2)));
            l_phase_maps_smooth(l,:) = circ_mean(ll); clear ll
            l=l+1;
           elseif size(firingMaps.rateMaps{condition},2) > 10 & ~isempty(fields{condition}{cell})
%             h_maps(h,:) = (squeeze(mean(firingMaps.rateMaps{condition}(cell,:,:),2)));
            for field = 1:length(fields{condition}{cell})
                hh = (squeeze((squeeze(binnedfiringMaps.phaseMaps{condition}(cell,:,:)))));
                for trial = 1:size(hh,1)
                   hh(trial,:) = circ_smoothTS(hh(trial,:),tau,'method','mean','exclude',0); 
                end
                h_rate_maps_smooth(h,:) = (squeeze(mean(firingMaps.rateMaps_box{condition}(cell,:,:),2)));
                h_phase_maps_smooth(h,:) = circ_mean(hh); clear hh
                COM(h) = fields{condition}{cell}{field}.COM;
                h=h+1;
            end
           end
            end
        end
    end
    
%     for ii=1:201
%     for jj=1:201
%     cc(ii,jj) = circ_corrcc(l_phase_maps_smooth(:,ii),l_phase_maps_smooth(:,jj));
%     hh(ii,jj) = circ_corrcc(h_phase_maps_smooth(:,ii),h_phase_maps_smooth(:,jj));
%     end
%     end
    subplot(4,2,1)
    imagesc(corr(h_rate_maps_smooth,'rows','complete'))
    caxis([0 1])
    subplot(4,2,2)
    imagesc(corr(l_rate_maps_smooth,'rows','complete'))
    caxis([0 1])
    subplot(4,2,3)
    imagesc(corr(h_phase_maps_smooth,'rows','complete'))
    caxis([0 1])
    subplot(4,2,4)
    imagesc(corr(l_phase_maps_smooth,'rows','complete'))
    caxis([0 1])
    subplot(4,2,5)
    imagesc(corr(cos(h_phase_maps_smooth),'rows','complete'))
    caxis([0 1])
    subplot(4,2,6)
    imagesc(corr(cos(l_phase_maps_smooth),'rows','complete'))
    caxis([0 1])
    subplot(4,2,7)
    imagesc(corr(sin(h_phase_maps_smooth),'rows','complete'))
    caxis([0 1])
%     imagesc(hh)
%     caxis([0 1])
    subplot(4,2,8)
    imagesc(corr(sin(l_phase_maps_smooth),'rows','complete'))
    caxis([0 1])
%     imagesc(cc)
%     caxis([0 1])
    pause(.01)
    
    end
%     cd D:\Dropbox\datasets\lsDataset
cd /home/david/datasets/lsDataset
end


% hp = circ_mean(h_phase_maps_smooth);
% lp = circ_mean(l_phase_maps_smooth);
% plot(hp)
% figure
% plot(lp,'g')







