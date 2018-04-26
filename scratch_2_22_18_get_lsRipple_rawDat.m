d = dir('*201*');
count=1;
count2=1;
for rec = 1:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    if exist([sessionInfo.FileName '.CA1Ripples.events.mat'])
    load([sessionInfo.FileName '.CA1Ripples.events.mat'])
%     bz_PlotRippleStats(ripples.maps,ripples.data,ripples.stats)
%     title([sessionInfo.FileName ])
%     pause
%     figure(1)
%     close
%     figure(2)
%     close
    
    %% get raw dat info
    warning off
    for ind = 1:size(ripples.timestamps,1)
        
data = bz_LoadBinary([sessionInfo.FileName '.dat'],'channels',ripples.rippleChan,'nChannels',length(sessionInfo.channels),'start', ripples.peaks(ind)-.15,'duration',.3);


% [a b] = min(data(2900:3100));
% offset = b-100;
% trough_shifted = circshift(data,-offset);
dat(count2,:) = data;%trough_shifted;
count2=count2+1;
% offsets(ind) = offset;
%         figure(100);
%         for i = 1:32
%             plot(squeeze(dat(ind,:,i))+i*1000,'k')
%             hold on
%         end
%         hold off
%         pause
    end
    
    %% get population firing rate
    spikes = bz_GetSpikes('noprompts',true);
    if ~isempty(spikes)
    for sp = 1:length(spikes.times)
        if strcmp(spikes.region{sp},'hpc') |strcmp(spikes.region{sp},'ca3') |strcmp(spikes.region{sp},'ca1')
            popRate(count,:) = zeros(501,1);
            for ind = 1:size(ripples.timestamps,1)
                ts = Restrict(spikes.times{sp},[ripples.peaks(ind)-.2 ripples.peaks(ind)+.2]); 
                idxNorm = ceil((ts-ripples.peaks(ind)+.2)*1250+.0000001);% - round(offsets(ind)./20000.*1250);
                if ~isempty(ts) & idxNorm > 0 & idxNorm <= 501
                    popRate(count,idxNorm) = popRate(count,idxNorm)+1;
                end
            end
            popRate(count,:) = popRate(count,:);
            popRate_z(count,:) = zscore(popRate(count,:));
            popRate_avg(count,:) = popRate(count,:)./size(ripples.timestamps,1);
            count = 1+count;
        end
    end
    end
    subplot(3,2,1)
    imagesc(popRate)
    subplot(3,2,2)
    plot(median(dat)*.195)
    subplot(3,2,3)
    plot((mean(popRate)))
    subplot(3,2,4)
    imagesc(popRate_avg)
    subplot(3,2,5)
    plot(mean(popRate_z))
    subplot(3,2,6)
    imagesc(popRate_z)
    caxis([0 4])
    pause(.1)
    end
   cd /home/david/datasets/ripples_LS 
save('/home/david/Dropbox/CA1Ripples_raw_traces.mat','-v7.3')
end