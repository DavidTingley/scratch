folders = dir('*201*');

for f=1:length(folders)
    cd(folders(f).name)
    if exist([folders(f).name '.spikes.cellinfo.mat']) & exist([folders(f).name '.CA1Ripples.events.mat'])
% 
%         f=1;
        sessionInfo = bz_getSessionInfo;
        load([sessionInfo.FileName '.behavior.mat'])
        spikes = bz_GetSpikes('noprompts',true);
        spkMat = bz_SpktToSpkmat(spikes,'binSize',.1,'overlap',4);
%        spkMat = bz_SpktToSpkmat(spikes,'binSize',.1,'overlap',1);
%        spkMat_large = bz_SpktToSpkmat(spikes,'binSize',.4,'overlap',4);
        hpc = find(strcmp(spikes.region ,'hpc') | strcmp(spikes.region ,'ca1'));
        latS = find(strcmp(spikes.region ,'ls'));
        reg(f) = 0;
        if isempty(hpc)
            hpc = find(strcmp(spikes.region ,'ca3'));
            reg(f) = 1;
        end

       
        %% carry on
        if ~isempty(hpc) & ~isempty(latS)
            for i=1:length(spikes.times)
                spkMat.dataZ(:,i) = zscore(spkMat.data(:,i));
            end
            hpcRates_z = nanmean(spkMat.dataZ(:,hpc)');
            hpcRates = nanmean(spkMat.data(:,hpc)');
            hpcPercentActive = nansum(spkMat.data(:,hpc)')./length(hpc);
            lsRates_z = nanmean(spkMat.dataZ(:,latS)');
            lsRates = nanmean(spkMat.data(:,latS)');
            
            hpcRates_smooth = fastrms(hpcRates,12);
            lsRates_smooth = fastrms(lsRates,12);
            hpcRates_smooth_z = fastrms(hpcRates_z,12)-fastrms(hpcRates_z,1200);
            lsRates_smooth_z = fastrms(lsRates_z,12);
          
            for i=1:95
                thresh_low = prctile(hpcRates_smooth,i-1);
                thresh_hi = prctile(hpcRates_smooth,i+5);
                idx = find(hpcRates_smooth>thresh_low & hpcRates_smooth<thresh_hi);
                rr(f,i) = nanmean(lsRates_smooth(idx));
                rz(f,i) = nanmean(lsRates_smooth_z(idx));
                
                thresh_low = prctile(hpcRates_smooth_z,i-1);
                thresh_hi = prctile(hpcRates_smooth_z,i+5);
                idx = find(hpcRates_smooth_z>thresh_low & hpcRates_smooth_z<thresh_hi);
                zz(f,i) = nanmean(lsRates_smooth_z(idx));
                zr(f,i) = nanmean(lsRates_smooth(idx));
                
            end
            nHPC_cells(f) = length(hpc);
            nLS_cells(f) = length(latS);
            
%             subplot(4,2,1)
%             imagesc(rr)
%             subplot(4,2,2)
%             imagesc(rz)
%             subplot(4,2,3)
%             imagesc(zr)
%             subplot(4,2,4)
            plot(nanmedian(zr(nHPC_cells>15,:)))
%             subplot(4,2,5)
%             plot(mean(zscore(ripHisto(nHPC_cells>15,:),[],2)))
%             hold on
%             plot(mean(zscore(thetaHisto(nHPC_cells>15,:),[],2)))
%             hold off            
%             subplot(4,2,6)
%             plot(mean(zscore(ripHisto_z(nHPC_cells>15,:),[],2)))
%             hold on
%             plot(mean(zscore(thetaHisto_z(nHPC_cells>15,:),[],2)))
%             hold off
%             subplot(4,2,7)
%             plot(mean((ripPercentHisto(nHPC_cells>15,:))))
%             hold on
%             plot(mean((thetaPercentHisto(nHPC_cells>15,:))))
%             hold off
            pause(.1)
        end
    
    end
    cd D:\datasets\lsDataset
end
% FileName = sessionInfo.FileName;
% clear lfp spikes spkMat* td_pow filt filt_lo phase ripples behavior sessionInfo
% save([FileName '.synchronyAnalysis.mat'],'-v7.3')
%% need to find all theta cycles and ripples and plot histograms over HPC synchrony percentiles
