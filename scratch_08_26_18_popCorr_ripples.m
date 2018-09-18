d = dir('*201*');
count=0;
ripWind = .1;
% [b a] = butter(4,[110/625 200/625],'bandpass');
[b a] = butter(4,[120/625 180/625],'bandpass');
    
%     
for rec = 1:length(d)
    cd(d(rec).name)
    if exist([d(rec).name '.behavior.mat']) & exist([d(rec).name '.firingMaps.cellinfo.mat']) & exist([d(rec).name '.CA1Ripples.events.mat'])
    load([d(rec).name '.behavior.mat'])
    load([d(rec).name '.firingMaps.cellinfo.mat'])
    load([d(rec).name '.CA1Ripples.events.mat'])
    spikes = bz_GetSpikes;
    sessionInfo = bz_getSessionInfo;
    
    hpc = find(strcmp(spikes.region,'hpc'));
    if ~isempty(hpc)
          
    for i=1:length(firingMaps.rateMaps)
        avgMap{i} = squeeze(mean(firingMaps.rateMaps{i},2));
        for j=1:length(spikes.times)
        avgMap_z{i}(j,:) = zscore(avgMap{i}(j,:));
        avgMap_zNorm{i}(j,:) = makeLength(avgMap_z{i}(j,:),ripWind * 400);
        avgMap_Norm{i}(j,:) = makeLength(avgMap{i}(j,:),ripWind * 400);
        end
    end

    for i=1:length(ripples.peaks)
        spkmat(i) = bz_SpktToSpkmat(spikes.times,'binSize',.05,'overlap',10,'win',[ripples.peaks(i)-ripWind ripples.peaks(i)+ripWind]);
    end
    for i=1:length(ripples.peaks)
        for j=1:length(spikes.times)
            spkmat(i).zscoredData(:,j) = zscore(spkmat(i).data(:,j));
        end
        b = spkmat(i).zscoredData(:,hpc);
%         b = b(1:80,:);
        for j=1:length(avgMap)
            co(j,i) = corr2(avgMap_zNorm{j}(hpc,:),b');
        end
    end
    start = find(min(behavior.events.trialIntervals(:))>ripples.peaks);
    start = start(end);
    stop = find(max(behavior.events.trialIntervals(:))<ripples.peaks);
    if isempty(stop)
        stop = length(ripples.peaks);
    else
        stop = stop(1);
    end
    subplot(2,2,1)
    imagesc(co)
    line([start start],[0 size(co,1)],'color','w')
    line([stop stop],[0 size(co,1)],'color','w')
    
    corrs{rec} = co; clear co;
    ints(rec,:) = [start stop];
    peaks{rec}= ripples.peaks;
    
        %% get LS lfp power too
    ls.lfp = bz_GetLFP(sessionInfo.ls,'noprompts', true);
    ls_power = (fastrms(bz_Filter(double(ls.lfp.data),'filter','butter','passband',[140 180],'order', 3),12));
    for event = 1:size(ripples.timestamps,1)
            start = round((ripples.peaks(event)-.02) * 1250);
            stop = round((ripples.peaks(event)+.02) * 1250);
            [ls_max(event) a] = max(abs(ls_power(start:stop)));
    end
    power{rec} = ls_max; clear ls_max;
    
    for i=1:size(corrs{rec},2)
        c(i) = max(corrs{rec}(:,i));
    end
    ls_power_hpc_replay(rec) = corr(c',power{rec}','rows','complete'); clear c
    ls_power_hpc_replay(ls_power_hpc_replay==0)=nan;
    subplot(2,2,4)
    histogram(ls_power_hpc_replay,40)
    pause(.1)
    end
    end
    cd /home/david/datasets/ripples_LS
    clear spkmat avgMap*
end

whos power peaks ints corrs
% save(['popCorr_ripples.mat'],'-v7.3')