folders = dir('*201*');

for f=1:length(folders)
    cd(folders(f).name)
    if exist([folders(f).name '.spikes.cellinfo.mat']) & exist([folders(f).name '.CA1Ripples.events.mat'])
        sessionInfo = bz_getSessionInfo;
        load([folders(f).name '.behavior.mat'])
        spikes = bz_GetSpikes('noprompts',true);
%         spkMat = bz_SpktToSpkmat(spikes,'binSize',.025,'overlap',1);
%         spkMat_large = bz_SpktToSpkmat(spikes,'binSize',.1,'overlap',4);
        spkMat = bz_SpktToSpkmat(spikes,'binSize',.1,'overlap',1);
        spkMat_large = bz_SpktToSpkmat(spikes,'binSize',.25,'overlap',5);
        hpc = find(strcmp(spikes.region ,'hpc'));
        reg(f) = 0;
        if isempty(hpc)
            hpc = find(strcmp(spikes.region ,'ca3'));
            reg(f) = 1;
        end
        if isempty(sessionInfo.ca1)
            lfpChan = sessionInfo.ca3;
        else
            lfpChan = sessionInfo.ca1;
        end
        %% get theta cycle times
        lfp = bz_GetLFP(lfpChan);
        [b a] = butter(3,[5/625 12/625],'bandpass');
        filt = filtfilt(b,a,double(lfp.data));
        [b a] = butter(3,[4/625],'low');
        filt_lo = filtfilt(b,a,double(lfp.data));
        td_pow = fastrms(filt,125)./fastrms(filt_lo,1250);
        phase = angle(hilbert(filt));
        
        [pks locs] = findpeaks(-phase);
        % cut by power thresh
        thresh = mean(td_pow)+std(td_pow)*3;
        idx = find(td_pow(locs)>thresh);
        cycles = lfp.timestamps(locs(idx));
        % cut by behavior
        id = find(InIntervals(cycles,[behavior.timestamps(1) behavior.timestamps(end)]));
        cycles = cycles(id);
        behavior.velocity(end+1)=behavior.velocity(end);
        for c =1:length(cycles)
            [a b]=min(abs(cycles(c)-behavior.timestamps));
            vels(c)=behavior.velocity(b);            
        end
        cycles(vels<20)=[]; clear vels
        
        thetaCycleCount = zeros(ceil(lfp.timestamps(end)*1000/25),1);
        thetaCycleCount(ceil(cycles*1000/25)) = 1;
        %% get ripple times
        load([sessionInfo.FileName '.CA1Ripples.events.mat'])
        rippleCount = zeros(ceil(lfp.timestamps(end)*1000/25),1);
        rippleCount(ceil(ripples.peaks*1000/25)) = 1;
        
        %% carry onn
        latS = find(strcmp(spikes.region ,'ls'));
        if ~isempty(hpc) & ~isempty(latS)
            for i=1:length(spikes.times)
                spkMat.dataZ(:,i) = zscore(spkMat.data(:,i));
            end
            hpcRates_z = nanmean(spkMat.dataZ(:,hpc)');
            hpcRates = nanmean(spkMat.data(:,hpc)');
            hpcPercentActive = nansum(spkMat_large.data(:,hpc)')./length(hpc);
            lsRates_z = nanmean(spkMat.dataZ(:,latS)');
            lsRates = nanmean(spkMat.data(:,latS)');
            
            hpcRates_smooth = fastrms(hpcRates,24);
            lsRates_smooth = fastrms(lsRates,24);
            hpcRates_smooth_z = fastrms(hpcRates_z,24);
            lsRates_smooth_z = fastrms(lsRates_z,24);
            for i=1:100
               idx = find(hpcPercentActive>i/100 & hpcPercentActive>(i+10)/100);
               ripPercentHisto(f,i) = nansum(rippleCount(idx));
               thetaPercentHisto(f,i) = nansum(thetaCycleCount(idx));
            end
            for i=1:95
                thresh_low = prctile(hpcRates_smooth,i-1);
                thresh_hi = prctile(hpcRates_smooth,i+5);
                idx = find(hpcRates_smooth>thresh_low & hpcRates_smooth<thresh_hi);
                rr(f,i) = nanmean(lsRates_smooth(idx));
                rz(f,i) = nanmean(lsRates_smooth_z(idx));
                
                ripHisto(f,i) = nansum(rippleCount(idx));
                thetaHisto(f,i) = nansum(thetaCycleCount(idx));
                
                thresh_low = prctile(hpcRates_smooth_z,i-1);
                thresh_hi = prctile(hpcRates_smooth_z,i+5);
                idx = find(hpcRates_smooth_z>thresh_low & hpcRates_smooth_z<thresh_hi);
                zz(f,i) = nanmean(lsRates_smooth_z(idx));
                zr(f,i) = nanmean(lsRates_smooth(idx));
                
                ripHisto_z(f,i) = nansum(rippleCount(idx));
                thetaHisto_z(f,i) = nansum(thetaCycleCount(idx));
            end
            nHPC_cells(f) = length(hpc);
            nLS_cells(f) = length(latS);
            
            subplot(4,2,1)
            imagesc(rr)
            subplot(4,2,2)
            imagesc(rz)
            subplot(4,2,3)
            imagesc(zr)
            subplot(4,2,4)
            plot(nanmean(zr(nHPC_cells>15,:)))
            subplot(4,2,5)
            plot(mean(zscore(ripHisto(nHPC_cells>15,:),[],2)))
            hold on
            plot(mean(zscore(thetaHisto(nHPC_cells>15,:),[],2)))
            hold off            
            subplot(4,2,6)
            plot(mean(zscore(ripHisto_z(nHPC_cells>15,:),[],2)))
            hold on
            plot(mean(zscore(thetaHisto_z(nHPC_cells>15,:),[],2)))
            hold off
            subplot(4,2,7)
            plot(mean(zscore(ripPercentHisto(nHPC_cells>15,:),[],2)))
            hold on
            plot(mean(zscore(thetaPercentHisto(nHPC_cells>15,:),[],2)))
            hold off
            pause(.1)
        end
    
    end
    cd D:\datasets\lsDataset
end

%% need to find all theta cycles and ripples and plot histograms over HPC synchrony percentiles
