d = dir('*201*');
count=0;
ls_rec =[];
hpc_rec= [];

animals = {'DT1','DT2','DT5','DT7','DT8','DT9'};
    
    
for rec = length(d):-1:1

    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo(pwd,'noprompts', true);
    if exist([sessionInfo.FileName '.LSRipples.events.mat']) & exist([sessionInfo.FileName '.CA1Ripples.events.mat'])
        ls = load([sessionInfo.FileName '.LSRipples.events.mat']);
        ca1 = load([sessionInfo.FileName '.CA1Ripples.events.mat']);
        [times groups] = spikes2sorted({ls.ripples.peaks,ca1.ripples.peaks});
        [ccg t] = CCG(times,groups,'binSize',.001);
        crossCorrs(rec,:) = ccg(:,1,2);
        if ~isfield(sessionInfo,'ls')
           sessionInfo.ls = ls.ripples.rippleChan;
           save([sessionInfo.FileName '.sessionInfo.mat'],'sessionInfo')
        end
        ls.lfp = bz_GetLFP(sessionInfo.ls,'noprompts', true);
        ca1.lfp = bz_GetLFP(sessionInfo.ca1,'noprompts', true);
%         [freqs,time,ls_power] = bz_WaveSpec(ls.lfp.data,[120 200],1,3,1/1250,'lin');
%         [freqs,time,hpc_power] = bz_WaveSpec(ca1.lfp.data,[120 200],1,3,1/1250,'lin');
%         ls_power = zscore(abs(ls_power));
%         hpc_power = zscore(abs(hpc_power));
ls_power = (fastrms(bz_Filter(double(ls.lfp.data),'filter','butter','passband',[140 180],'order', 3),12));
hpc_power = (fastrms(bz_Filter(double(ca1.lfp.data),'filter','butter','passband',[140 180],'order', 3),12));
        
figure(1)
idx = find(strcmp(animals,sessionInfo.animal));
        subplot(6,3,idx);
        for event = 1:size(ls.ripples.timestamps,1)
            start = round((ls.ripples.peaks(event)-.02) * 1250);
            stop = round((ls.ripples.peaks(event)+.02) * 1250);

            [ls_max a] = max(abs(ls_power(start:stop)));
            [hpc_max a] = max(abs(hpc_power(start:stop)));
            scatter(ls_max,hpc_max,2,'.r')
            ls_rec=[ls_rec;ls_max,hpc_max];
            hold on
        end
        for event = 1:size(ca1.ripples.timestamps,1)
            start = round((ca1.ripples.peaks(event)-.02) * 1250);
            stop = round((ca1.ripples.peaks(event)+.02) * 1250);

            [ls_max a] = max(abs(ls_power(start:stop)));
            [hpc_max a] = max(abs(hpc_power(start:stop)));
            scatter(ls_max,hpc_max,2,'.k')
            hpc_rec=[hpc_rec;ls_max,hpc_max];
            hold on
        end

        subplot(6,3,idx+6);
        pts = linspace(0, max([hpc_rec(:)]), 101);
        N = histcounts2(hpc_rec(:,1), hpc_rec(:,2), pts, pts);
        imagesc(pts,pts,log(N))
        title('hpc detection')
        subplot(6,3,idx+12);
        N = histcounts2(ls_rec(:,1), ls_rec(:,2), pts, pts);
        imagesc(pts,pts,log(N))
        title('ls detection')
        
        pause(1)
%% get population firing rate
figure(2)
    spikes = bz_GetSpikes('noprompts',true);
    if ~isempty(spikes)
    for sp = 1:length(spikes.times)
        if strcmp(spikes.region{sp},'ls') %|strcmp(spikes.region{sp},'ca3') |strcmp(spikes.region{sp},'ca1')
            popRate(count,:) = zeros(501,1);
            for ind = 1:size(ripples.timestamps,1)
                ts = Restrict(spikes.times{sp},[ca1.ripples.peaks(ind)-.2 ca1.ripples.peaks(ind)+.2]); 
                idxNorm = ceil((ts-ca1.ripples.peaks(ind)+.2)*1250+.0000001);% - round(offsets(ind)./20000.*1250);
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
    plot(squeeze(median(median(dat,2)))*.195)
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
%    cd /home/david/datasets/ripples_LS 
cd E:\datasets\ripples_LS
end