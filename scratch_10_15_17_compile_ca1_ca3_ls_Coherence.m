d = dir('*201*');
coh_ca3 = zeros(1,1000);
coh_ca1 = zeros(1,1000);
c1=1;
c3=1;
[b a] = butter(4,[4/625 12/625],'bandpass');

for ii=1:61
    cd(d(ii).name)


    sessionInfo = bz_getSessionInfo;
%     if ~isempty(sessionInfo.ca3) & ~isempty(sessionInfo.ca1)
    load([sessionInfo.FileName '.behavior.mat']);
    if strcmp(behavior.description,'wheel alternation')
    ca1 = [];
    ca3 = [];
    ls = [];
    if ~isempty(sessionInfo.ca1)
        ca1 = bz_GetLFP(sessionInfo.ca1);
        ca1.data = FiltFiltM(b,a,double(ca1.data));
    end
    if ~isempty(sessionInfo.ca3)
        ca3 = bz_GetLFP(sessionInfo.ca3);
        ca3.data = FiltFiltM(b,a,double(ca3.data));
    end
    if ~isempty(sessionInfo.ls)
        ls = bz_GetLFP(sessionInfo.ls);
        ls.data = FiltFiltM(b,a,double(ls.data));
    end
%     save([sessionInfo.FileName '.lfp.mat'],'ls','ca1','ca3')

    if ~isempty(sessionInfo.ca1)
        for i=1:length(behavior.events.trialConditions)
            start = round(behavior.events.trialIntervals((i),1)*1250);
            stop = round(behavior.events.trialIntervals((i),2)*1250);
            [coherogram phase t f] = bz_MTCoherogram(double(ls.data(start:stop,1)),double(ca1.data(start:stop,1)),'range',[4 12],'window',1,'step',1./1250);
            coh_ca1(c1,:) = makeLength(mean(coherogram),1000);
            c1=c1+1;
        end
    end
    if ~isempty(sessionInfo.ca3)
        for i=1:length(behavior.events.trialConditions)
            start = round(behavior.events.trialIntervals((i),1)*1250);
            stop = round(behavior.events.trialIntervals((i),2)*1250);
            [coherogram phase t f] = bz_MTCoherogram(double(ls.data(start:stop,1)),double(ca3.data(start:stop,1)),'range',[4 12],'window',1,'step',1./1250);
            coh_ca3(c3,:) = makeLength(mean(coherogram),1000);
            c3=c3+1;
        end
    end
    subplot(2,2,1)
    plot(nanmean(coh_ca3),'r')
    hold on
    plot(nanmean(coh_ca1),'g')
    hold off
    subplot(2,2,2)
    plot(nanmedian(coh_ca3),'r')
    hold on
    plot(nanmedian(coh_ca1),'g')
    hold off
    subplot(2,2,3)
    imagesc(coh_ca3)
    subplot(2,2,4)
    imagesc(coh_ca1)
    pause(.01)
%     end
    end
    clear ls ca1 ca3
    cd /home/david/datasets/lsDataset
end