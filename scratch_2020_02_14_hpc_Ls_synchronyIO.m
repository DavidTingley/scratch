% folders = dir('*201*');
% 
% for f=1:length(folders)
%     cd(folders(f).name)
%     if exist([folders(f).name '.spikes.cellinfo.mat']) & exist([folders(f).name '.CA1Ripples.events.mat'])

        f=1;
        sessionInfo = bz_getSessionInfo;
        load([sessionInfo.FileName '.behavior.mat'])
        spikes = bz_GetSpikes('noprompts',true);
%         spkMat = bz_SpktToSpkmat(spikes,'binSize',.1,'overlap',4);
        spkMat = bz_SpktToSpkmat(spikes,'binSize',.25,'overlap',10);
%        spkMat_large = bz_SpktToSpkmat(spikes,'binSize',.4,'overlap',4);
        hpc = find(strcmp(spikes.region ,'hpc') | strcmp(spikes.region ,'ca1'));
        latS = find(strcmp(spikes.region ,'ls'));
        reg(f) = 0;
        if isempty(hpc)
            hpc = find(strcmp(spikes.region ,'ca3'));
            reg(f) = 1;
        end
        if ~isempty(sessionInfo.ca3)
            lfpChan = sessionInfo.ca3;
        else
            lfpChan = sessionInfo.ca1;
        end
        %% get theta cycle times
        lfp = bz_GetLFP(lfpChan);
        if exist([sessionInfo.FileName '.CA1Ripples.events.mat'])
            [b a] = butter(4,[5/625 12/625],'bandpass');
            filt = filtfilt(b,a,double(lfp.data));
            [b a] = butter(4,[4/625],'low');
            filt_lo = filtfilt(b,a,double(lfp.data));
            td_pow = fastrms(filt,125)./fastrms(filt_lo,1250);
            phase = angle(hilbert(filt));

            [pks locs] = findpeaks(-phase,'MinPeakDistance',100);
            % cut by power thresh
            powThresh(f) = mean(td_pow)+std(td_pow);
            idx = find(td_pow(locs)>powThresh(f));
            cycles = lfp.timestamps(locs(idx));
            % cut by behavior
            id = find(InIntervals(cycles,[behavior.events.trialIntervals]));
            cycles = cycles(id);
            behavior.velocity(end+1)=behavior.velocity(end);

%             vels=behavior.velocity(ceil((cycles-behavior.timestamps(1))*round(behavior.samplingRate)));
            parfor c =1:length(cycles)
                [a b]=min(abs(cycles(c)-behavior.timestamps));
                vels(c)=behavior.velocity(b);            
            end
            cycles(vels<20)=[]; 
            velocities{f} = vels; clear vels

            powers{f} = td_pow(round(cycles*1250));   
%             
            %% get theta cycle times
            thetaCycleCount = zeros(ceil(lfp.timestamps(end)*1000/25),1);
            thetaCycleCount(ceil(cycles*1000/25)) = 1;
            %% get ripple times
            load([sessionInfo.FileName '.CA1Ripples.events.mat'])
            rippleCount = zeros(ceil(lfp.timestamps(end)*1000/25),1);
            rippleCount(ceil(ripples.peaks*1000/25)) = 1;
            %% spikes per event analysis
            thresh = .06;
            spkIndices = spikes.spindices(ismember(spikes.spindices(:,2),hpc));
            spkUids = spikes.spindices(ismember(spikes.spindices(:,2),hpc),2);
            spkIndices_ls = spikes.spindices(ismember(spikes.spindices(:,2),latS));
            tempPETH = zeros(length(cycles),120);
            parfor cyc = 1:length(cycles)
                id = FindInInterval(spkIndices,[cycles(cyc)-thresh cycles(cyc)+thresh]);
                if ~isempty(id)
                    temp(cyc) = diff(id)+1;
                    s = round((spkIndices(id(1):id(2))-cycles(cyc)+.06)*1000);
                    tempPETH(cyc,:) = hist(s,1:120);
                    tempCells(cyc) = length(unique(spkUids(id(1):id(2))));
                else
                    temp(cyc) = 0;
                    tempPeth(cyc,:) = zeros(1,120);
                    tempCells(cyc) = 0;
                end
                
                if ~isempty(spkIndices_ls) & ~isempty(FindInInterval(spkIndices_ls,[cycles(cyc)-thresh cycles(cyc)+thresh]))
                    temp2(cyc) = diff(FindInInterval(spkIndices_ls,[cycles(cyc)-thresh cycles(cyc)+thresh]))+1;
                end
            end
            cellsPerThetaCyc{f} = tempCells; clear tempCells
            thetaCycPETH{f} = tempPETH; clear tempPETH
            spksPerThetaCyc{f} = temp; clear temp
            if exist('temp2')
            spksPerThetaCyc_ls{f} = temp2; clear temp2
            end
            tempPETH = zeros(length(ripples.peaks),120);
            parfor rip =1:length(ripples.peaks)
                id = FindInInterval(spkIndices,[ripples.peaks(rip)-thresh ripples.peaks(rip)+thresh]);
                if ~isempty(id)
                    temp(rip) = diff(id)+1;
                    s = round((spkIndices(id(1):id(2))-ripples.peaks(rip)+.06)*1000);
                    tempPETH(rip,:) = hist(s,1:120);
                    tempRip(rip) = length(unique(spkUids(id(1):id(2))));
                else
                    temp(rip) = 0;
                    tempPeth(rip,:) = zeros(1,120);
                    tempRip(rip) = 0;
                end
                if ~isempty(spkIndices_ls) & ~isempty(FindInInterval(spkIndices_ls,[ripples.peaks(rip)-thresh ripples.peaks(rip)+thresh]))
                    temp2(rip) = diff(FindInInterval(spkIndices_ls,[ripples.peaks(rip)-thresh ripples.peaks(rip)+thresh]))+1;
                end
            end
            cellsPerRipple{f} = tempRip; clear tempRip
            ripPETH{f} = tempPETH; clear tempPETH
            spksPerRipple{f} = temp; clear temp
            if exist('temp2')
            spksPerRipple_ls{f} = temp2; clear temp2
            end
        else
            rippleCount = nan(ceil(lfp.timestamps(end)*1000/25),1);
            thetaCycleCount = nan(ceil(lfp.timestamps(end)*1000/25),1);
            ripPETH{f} = nan(1,120);
            thetaCycPETH{f} = nan(1,120);
            spksPerRipple_ls{f}=[];
            spksPerRipple{f}=[];
            spksPerThetaCyc_ls{f}=[];
            spksPerThetaCyc{f}=[];
        end
        %% carry on
        if ~isempty(hpc) & ~isempty(latS)
            for i=1:length(spikes.times)
                spkMat.dataZ(:,i) = zscore(spkMat.data(:,i));
            end
            hpcCounts = nansum(spkMat.data(:,hpc)'>0);
            hpcRates_z = nanmean(spkMat.dataZ(:,hpc)');
            hpcRates = nanmean(spkMat.data(:,hpc)');
            hpcPercentActive = nansum(spkMat.data(:,hpc)')./length(hpc);
            lsRates_z = nanmean(spkMat.dataZ(:,latS)');
            lsRates = nanmean(spkMat.data(:,latS)');
            
            hpcCounts_smooth = fastrms(hpcCounts,24);
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
                
                thresh_low = prctile(hpcCounts_smooth,i-1);
                thresh_hi = prctile(hpcCounts_smooth,i+5);
                idx = find(hpcCounts_smooth>thresh_low & hpcCounts_smooth<thresh_hi);
                ripHisto_count(f,i) = nansum(rippleCount(idx));
                thetaHisto_count(f,i) = nansum(thetaCycleCount(idx));
                count_z(f,i) = nanmean(lsRates_smooth_z(idx));
                count_r(f,i) = nanmean(lsRates_smooth(idx));
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
%             plot(nanmean(zr(nHPC_cells>15,:)))
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
%             pause(.1)
        end
    
%     end
%     cd D:\datasets\lsDataset
% end
FileName = sessionInfo.FileName;
clear lfp spikes spkMat* td_pow filt filt_lo phase ripples behavior sessionInfo
save([FileName '.synchronyAnalysis.mat'],'-v7.3')
%% need to find all theta cycles and ripples and plot histograms over HPC synchrony percentiles
