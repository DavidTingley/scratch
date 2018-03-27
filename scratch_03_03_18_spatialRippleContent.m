d = dir('*201*');
count=0;
[b a] = butter(4,[110/625 200/625],'bandpass');
% [b a] = butter(4,[120/625 180/625],'bandpass');
    
%     
% for rec = 1:length(d)
%     cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    if  exist([sessionInfo.FileName '.CA1Ripples.events.mat']) & exist([sessionInfo.FileName '.olypherInfo_w_disc.cellinfo.mat']) & exist([sessionInfo.FileName '.olypherInfo_w_disc.cellinfo.mat'])
        
        load([sessionInfo.FileName '.olypherInfo_w_disc.cellinfo.mat'],'olypherInfo')
        load([sessionInfo.FileName '.rewardModulation.cellinfo.mat'],'rewardModulation')
        ca1 = load([sessionInfo.FileName '.CA1Ripples.events.mat']);
        spikes = bz_GetSpikes;
        
        %% get LFP power for all events..        
        ls.lfp = bz_GetLFP(sessionInfo.ls);
        ca1.lfp = bz_GetLFP(sessionInfo.ca1);
        ls_power = (fastrms(FiltFiltM(b,a,double(ls.lfp.data)),15));
        ls_power_z = zscore(ls_power);
        hpc_power = (fastrms(FiltFiltM(b,a,double(ca1.lfp.data)),15));
        ls_rec =[];
        hpc_rec= [];
        for event = 1:size(ca1.ripples.timestamps,1)
            start = round((ca1.ripples.timestamps(event,1)-.015) * 1250); % used to be 20 ms
            stop = round((ca1.ripples.timestamps(event,1)+.015) * 1250);

            [ls_max blah] = max(abs(ls_power(start:stop)));
            [ls_max_z blah_z] = max(abs(ls_power_z(start:stop)));
            [hpc_max blah] = max(abs(hpc_power(start:stop)));
%             scatter(ls_max,hpc_max,'.k')
            hpc_rec=[hpc_rec;ls_max,hpc_max];
            ls_rec = [ls_rec;ls_max_z];
%             hold on
        end
        
        %% get that content
        spatialContent = zeros(length(spikes.times),length(ca1.ripples.peaks));
        rewardContent = zeros(length(spikes.times),length(ca1.ripples.peaks));
        nSpikes = zeros(length(spikes.times),length(ca1.ripples.peaks));
        spatialContent_part = zeros(length(spikes.times),length(ca1.ripples.peaks));
        rewardContent_part = zeros(length(spikes.times),length(ca1.ripples.peaks));
        
        for spk = 1:length(spikes.times)
            if strcmp(spikes.region{spk},'hpc') | strcmp(spikes.region{spk},'ca3') | strcmp(spikes.region{spk},'ca1') 
            rows = find(olypherInfo.results{spk}.discBins==2);
            cols = find(olypherInfo.results{spk}.smoothing==20);
            cols = intersect(rows,cols);

            meanPeakRate(spk) = nanmean(olypherInfo.results{spk}.ratePeakInfo(cols));
            for ind = 1:length(ca1.ripples.peaks)
                ripSpks = Restrict(spikes.times{spk},[ca1.ripples.timestamps(ind,1)-.01 ca1.ripples.timestamps(ind,2)+.01]);
                spatialContent(spk,ind) = meanPeakRate(spk).*length(ripSpks);
                rewardContent(spk,ind) = rewardModulation.rewardGain(spk).*length(ripSpks);
                nSpikes(spk,ind) = length(ripSpks);
                if length(ripSpks) > 0
                spatialContent_part(spk,ind) = meanPeakRate(spk);%.*length(ripSpks);
                rewardContent_part(spk,ind) = rewardModulation.rewardGain(spk);%.*length(ripSpks);
                end
            end
            count = 1+count;      
           
            else 
                meanPeakRate = [];
            end
%             spk
%             subplot(3,2,1)
%             scatter(mean(spatialContent),mean(rewardContent),'.r')
%             subplot(3,2,2)
%             scatter(mean(spatialContent_part),mean(rewardContent_part),'.r')
%             subplot(3,2,3)
%             scatter(mean(spatialContent),hpc_rec(:,1),'.r')
%             subplot(3,2,4)
%             scatter(mean(rewardContent),hpc_rec(:,1),'.r')
%             subplot(3,2,5)
%             scatter(mean(spatialContent_part),hpc_rec(:,1),'.r')
%             subplot(3,2,6)
%             scatter(mean(rewardContent_part),hpc_rec(:,1),'.r')
%             pause(.1)
        end
    content.region{rec} = spikes.region{end};
    content.nCells(rec) = size(spatialContent,1);
    content.spatialContent{rec} = spatialContent;
    content.rewardContent{rec} = rewardContent;
    content.spatialContent_binary{rec} = spatialContent_part;
    content.rewardContent_binary{rec} = rewardContent_part;
    content.nSpikes{rec} = nSpikes;
    content.meanPeakRate{rec} = meanPeakRate;
    content.rewardGain{rec} = rewardModulation.rewardGain;
    content.hpc_power{rec} = hpc_rec(:,2);
    content.ls_power{rec} = hpc_rec(:,1);
    content.ls_power_z{rec} = ls_rec(:,1);
    [aa bb] = corr(content.ls_power{rec},nanmean(content.spatialContent{rec})')
    [aa bb] = corr(content.ls_power{rec},nanmean(content.rewardContent{rec})')
    rec
    clear meanPeakRate rewardModulation hpc_rec ls_rec meanPeakRate nSpikes *Content*
    end
    save([sessionInfo.FileName '.rippleContent.mat'])
%    cd /home/david/datasets/ripples_LS 
%    save('/home/david/Dropbox/hpc_ripple_content_180324.mat','-v7.3')
% end