d = dir('*201*');
count=0;
    ls_rec =[];
    hpc_rec= [];
for rec = length(d):-1:1

    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    if exist([sessionInfo.FileName '.LSRipples.events.mat']) & exist([sessionInfo.FileName '.CA1Ripples.events.mat'])
        ls = load([sessionInfo.FileName '.LSRipples.events.mat']);
        ca1 = load([sessionInfo.FileName '.CA1Ripples.events.mat']);
        count = 1+count;
        ls.lfp = bz_GetLFP(sessionInfo.ls);
        ca1.lfp = bz_GetLFP(sessionInfo.ca1);
%         [freqs,time,ls_power] = bz_WaveSpec(ls.lfp.data,[120 200],1,3,1/1250,'lin');
%         [freqs,time,hpc_power] = bz_WaveSpec(ca1.lfp.data,[120 200],1,3,1/1250,'lin');
%         ls_power = zscore(abs(ls_power));
%         hpc_power = zscore(abs(hpc_power));
ls_power = (fastrms(bz_Filter(double(ls.lfp.data),'filter','butter','passband',[120 180],'order', 3),20));
hpc_power = (fastrms(bz_Filter(double(ca1.lfp.data),'filter','butter','passband',[120 180],'order', 3),20));
        
%         figure(rec)
        subplot(2,2,1);
        for event = 1:size(ls.ripples.timestamps,1)
            start = round((ls.ripples.peaks(event)-.02) * 1250);
            stop = round((ls.ripples.peaks(event)+.02) * 1250);

            [ls_max a] = max(abs(ls_power(start:stop)));
            [hpc_max a] = max(abs(hpc_power(start:stop)));
            scatter(ls_max,hpc_max,'.r')
            ls_rec=[ls_rec;ls_max,hpc_max];
            hold on
        end
        for event = 1:size(ca1.ripples.timestamps,1)
            start = round((ca1.ripples.peaks(event)-.02) * 1250);
            stop = round((ca1.ripples.peaks(event)+.02) * 1250);

            [ls_max a] = max(abs(ls_power(start:stop)));
            [hpc_max a] = max(abs(hpc_power(start:stop)));
            scatter(ls_max,hpc_max,'.k')
            hpc_rec=[hpc_rec;ls_max,hpc_max];
            hold on
        end
%         axis([0 40 0 40])
        
        subplot(2,2,2);
        pts = linspace(0, max([hpc_rec(:)]), 101);
        N = histcounts2(hpc_rec(:,1), hpc_rec(:,2), pts, pts);
        imagesc(pts,pts,log(N))
        title('hpc detection')
        subplot(2,2,3);
        N = histcounts2(ls_rec(:,1), ls_rec(:,2), pts, pts);
        imagesc(pts,pts,log(N))
        title('ls detection')
        
        pause(1)
%                 ls.dat = bz_LoadBinary([sessionInfo.FileName '.dat'],'nChannels',length(sessionInfo.channels),'channels',sessionInfo.ls);
%         ca1.dat = bz_LoadBinary([sessionInfo.FileName '.dat'],'nChannels',length(sessionInfo.channels),'channels',sessionInfo.ca1);
%         
%         for event = 1:size(ls.ripples.timestamps,1)
%             start = round((ls.ripples.peaks(event)-.02) * 20000);
%             stop = round((ls.ripples.peaks(event)+.02) * 20000);
%             [freqs,time,ls_power] = bz_WaveSpec(ls.dat.data(start:stop),[1 200],200,4,1/20000,'lin');
%             [freqs,time,hpc_power] = bz_WaveSpec(ca1.dat.data(start:stop),[1 200],200,4,1/20000,'lin');
%             [a ls_max] = max(abs(ls_power));
%             [a hpc_max] = max(abs(hpc_power));
%             scatter(ls_max,hpc_max,'.r')
%             hold on
%         end
        
%         [times groups] = spikes2sorted({ls.ripples.peaks,ca1.ripples.peaks});
%         [ccg t] = CCG(times,groups,'binSize',.001);
    end
   cd /home/david/datasets/ripples_LS 
end