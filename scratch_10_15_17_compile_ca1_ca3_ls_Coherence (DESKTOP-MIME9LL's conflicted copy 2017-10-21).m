d = dir('*201*');
coh_ca3 = zeros(1,1000);
coh_ca1 = zeros(1,1000);
coh_ca3_phase = zeros(1,1000);
coh_ca1_phase = zeros(1,1000);
c1=1;
c3=1;
[b a] = butter(3,[4/625 12/625],'bandpass');
for ii=1:length(d)
    cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
%     if ~isempty(sessionInfo.ca3) & ~isempty(sessionInfo.ca1)
    load([sessionInfo.FileName '.behavior.mat']);
    if exist([sessionInfo.FileName '.lfp.mat']) %strcmp(behavior.description,'wheel alternation')
%     ca1 = [];
%     ca3 = [];
%     ls = [];
%     if ~isempty(sessionInfo.ca1)
%         ca1 = bz_GetLFP(sessionInfo.ca1);
%         ca1.data = FiltFiltM(b,a,double(ca1.data));
%     end
%     if ~isempty(sessionInfo.ca3)
%         ca3 = bz_GetLFP(sessionInfo.ca3);
%         ca3.data = FiltFiltM(b,a,double(ca3.data));
%     end
%     if ~isempty(sessionInfo.ls)
%         ls = bz_GetLFP(sessionInfo.ls);
%         ls.data = FiltFiltM(b,a,double(ls.data));
%     end
    load([sessionInfo.FileName '.lfp.mat'],'ls','ca1','ca3')

    if ~isempty(sessionInfo.ca1)
        for i=1:length(behavior.events.trialConditions)
            if strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'linear')
                ind=1;
            elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'central')
                ind=2;
            elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'wheel')
                ind=3;
            end
            start = round(behavior.events.trialIntervals((i),1)*1250);
            stop = round(behavior.events.trialIntervals((i),2)*1250);
            if stop-start > 1000
            [coherogram phase t f] = bz_MTCoherogram(double(ls.data(start:stop,1)),double(ca1.data(start:stop,1)),'range',[4 12],'window',1,'step',1./1250);
            coh_ca1(c1,:) = makeLength(mean(coherogram),1000);
            [coherogram phase t f] = bz_MTCoherogram(angle(hilbert(FiltFiltM(b,a,double(ls.data(start:stop,1))))),angle(hilbert(FiltFiltM(b,a,double(ca1.data(start:stop,1))))),'range',[4 12],'window',1,'step',1./1250);
            coh_ca1_phase(c1,:) = makeLength(mean(coherogram),1000);
            ca1_ind(c1) = ind;
            c1=c1+1;
            end
        end
    end
    if ~isempty(sessionInfo.ca3)
        for i=1:length(behavior.events.trialConditions)
            if strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'linear')
                ind=1;
            elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'central')
                ind=2;
            elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'wheel')
                ind=3;
            end
            start = round(behavior.events.trialIntervals((i),1)*1250);
            stop = round(behavior.events.trialIntervals((i),2)*1250);
            if stop-start > 1000
            [coherogram phase t f] = bz_MTCoherogram(double(ls.data(start:stop,1)),double(ca3.data(start:stop,1)),'range',[4 12],'window',1,'step',1./1250);
            coh_ca3(c3,:) = makeLength(mean(coherogram),1000);
            [coherogram phase t f] = bz_MTCoherogram(angle(hilbert(FiltFiltM(b,a,double(ls.data(start:stop,1))))),angle(hilbert(FiltFiltM(b,a,double(ca3.data(start:stop,1))))),'range',[4 12],'window',1,'step',1./1250);
            coh_ca3_phase(c3,:) = makeLength(mean(coherogram),1000);
            ca3_ind(c3) = ind;
            c3=c3+1;
            end
        end
    end
    for i=1:3
        figure(i)
    subplot(2,2,1)
    plot(nanmean(coh_ca3(ca3_ind==i,:)),'r')
    hold on
    plot(nanmean(coh_ca1(ca1_ind==i,:)),'g')
    hold off
    subplot(2,2,2)
    plot(nanmean(coh_ca3_phase(ca3_ind==i,:)),'r')
    hold on
    plot(nanmean(coh_ca1_phase(ca1_ind==i,:)),'g')
    hold off
    subplot(2,2,3)
    imagesc(coh_ca3(ca3_ind==i,:))
    subplot(2,2,4)
    imagesc(coh_ca1(ca1_ind==i,:))
    end
    pause(.01)
%     end
    end
    clear ls ca1 ca3
    cd('D:\Dropbox\datasets\lsDataset')
%     cd /home/david/datasets/lsDataset
end