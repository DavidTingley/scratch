d = dir('*201*');
coh_ca3 = zeros(1,1000);
coh_ca1 = zeros(1,1000);
coh_ca3_phase = zeros(1,1000);
coh_ca1_phase = zeros(1,1000);
coh_ca3_z = zeros(1,1000);
coh_ca1_z = zeros(1,1000);
coh_ca3_phase_z = zeros(1,1000);
coh_ca1_phase_z = zeros(1,1000);
ca1_ind = [];
ca3_ind = [];
coh_ca3_theta_gamma_z = zeros(1,1000);
coh_ca1_theta_gamma_z = zeros(1,1000);
coh_ca1_ca3 = zeros(1,1000);
coh_ca1_ca3_phase = zeros(1,1000);
coh_ca1_ca3_z = zeros(1,1000);
coh_ca1_ca3_phase_z = zeros(1,1000);
ca1_ca3_ind = [];
ca1_recording=[];ca3_recording=[];
c1=1;
c2=1;
c3=1;
[b a] = butter(3,[4/625 12/625],'bandpass');

[bbb aaa] = butter(3,[40/625 200/625],'bandpass');
for ii=1:length(d)
    cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    disp(['working on: ' sessionInfo.FileName])
    if exist([sessionInfo.FileName '.lfp'])
    
%     if ~isempty(sessionInfo.ca3) & ~isempty(sessionInfo.ca1)
    load([sessionInfo.FileName '.behavior.mat']);
%     if strcmp(behavior.description,'wheel alternation')
    ca1 = [];
    ca3 = [];
    ls = [];
    if ~isempty(sessionInfo.ca1)
        ca1 = bz_GetLFP(sessionInfo.ca1);
%         ca1.data = FiltFiltM(b,a,double(ca1.data));
    end
    if ~isempty(sessionInfo.ca3)
        ca3 = bz_GetLFP(sessionInfo.ca3);
%         ca3.data = FiltFiltM(b,a,double(ca3.data));
    end
    if ~isempty(sessionInfo.ls)
        ls = bz_GetLFP(sessionInfo.ls);
%         ls.data = FiltFiltM(b,a,double(ls.data));
    end
%     load([sessionInfo.FileName '.lfp.mat'],'ls','ca1','ca3')

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
            [coherogram phase t f] = bz_MTCoherogram(double(ls.data(start:stop,1)),double(ca1.data(start:stop,1)),'range',[4 12],'window',.15,'step',.05);
            coh_ca1(c1,:) = makeLength(mean(coherogram),1000);
            [coherogram phase t f] = bz_MTCoherogram(angle(hilbert(FiltFiltM(b,a,double(ls.data(start:stop,1))))),angle(hilbert(FiltFiltM(b,a,double(ca1.data(start:stop,1))))),'range',[4 12],'window',.15,'step',.05);
            coh_ca1_phase(c1,:) = makeLength(mean(coherogram),1000);
            coh_ca1_z(c1,:) = zscore(coh_ca1(c1,:));
            coh_ca1_phase_z(c1,:) = zscore(coh_ca1_phase(c1,:));
            
%             t = fastrms(FiltFiltM(bbb,aaa,double(ls.data(start:stop,1))),60);
%             tt = angle(hilbert(FiltFiltM(b,a,double(ca1.data(start:stop,1)))));
%             for temp = 1:length(t)
%                 if temp > length(t) - 300
% %                     c(temp) = circ_corrcl(tt(temp-150:end),t(temp-150:end));
%                 elseif temp < 301
% %                     c(temp) = circ_corrcl(tt(1:temp+150),t(1:temp+150));
%                 else
%                 	c(temp) = circ_corrcl(tt(temp-300:temp+300),t(temp-300:temp+300));
%                 end
%             end
            [coherogram phase t f] = bz_MTCoherogram(fastrms(FiltFiltM(bbb,aaa,double(ls.data(start:stop,1))),20),angle(hilbert(FiltFiltM(b,a,double(ca1.data(start:stop,1))))),'range',[4 12],'window',.3,'step',.05);
            coh_ca1_theta_gamma(c1,:) = makeLength(mean(coherogram),1000); %makeLength(c(301:end),1000); clear c
            coh_ca1_theta_gamma_z(c1,:) = zscore(coh_ca1_theta_gamma(c1,:));
            ca1_animal(c1) = sum(double(sessionInfo.animal));
            ca1_recording(c1) = ii;
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
            [coherogram phase t f] = bz_MTCoherogram(double(ls.data(start:stop,1)),double(ca3.data(start:stop,1)),'range',[4 12],'window',.15,'step',.05);
            coh_ca3(c3,:) = makeLength(mean(coherogram),1000);
            [coherogram phase t f] = bz_MTCoherogram(angle(hilbert(FiltFiltM(b,a,double(ls.data(start:stop,1))))),angle(hilbert(FiltFiltM(b,a,double(ca3.data(start:stop,1))))),'range',[4 12],'window',.15,'step',.05);
            coh_ca3_phase(c3,:) = makeLength(mean(coherogram),1000);
            coh_ca3_z(c3,:) = zscore(coh_ca3(c3,:));
            coh_ca3_phase_z(c3,:) = zscore(coh_ca3_phase(c3,:));
            
%             t = fastrms(FiltFiltM(bbb,aaa,double(ls.data(start:stop,1))),60);
%             tt = angle(hilbert(FiltFiltM(b,a,double(ca3.data(start:stop,1)))));
%             for temp = 1:length(t)
%                 if temp > length(t) - 300
% %                     c(temp) = circ_corrcl(tt(temp-150:end),t(temp-150:end));
%                 elseif temp < 301
% %                     c(temp) = circ_corrcl(tt(1:temp+150),t(1:temp+150));
%                 else
%                 	c(temp) = circ_corrcl(tt(temp-300:temp+300),t(temp-300:temp+300));
%                 end
%             end
            [coherogram phase t f] = bz_MTCoherogram(fastrms(FiltFiltM(bbb,aaa,double(ls.data(start:stop,1))),20),angle(hilbert(FiltFiltM(b,a,double(ca3.data(start:stop,1))))),'range',[4 12],'window',.3,'step',.05);
            coh_ca3_theta_gamma(c3,:) =  makeLength(mean(coherogram),1000); %makeLength(c(301:end),1000); clear c
            coh_ca3_theta_gamma_z(c3,:) = zscore(coh_ca3_theta_gamma(c3,:));
            ca3_recording(c3) = ii;
            ca3_animal(c3) = sum(double(sessionInfo.animal));
            ca3_ind(c3) = ind;
            c3=c3+1;
            end
        end
    end
    if ~isempty(sessionInfo.ca3) & ~isempty(sessionInfo.ca1)
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
            [coherogram phase t f] = bz_MTCoherogram(double(ca1.data(start:stop,1)),double(ca3.data(start:stop,1)),'range',[4 12],'window',.15,'step',.05);
            coh_ca1_ca3(c2,:) = makeLength(mean(coherogram),1000);
            [coherogram phase t f] = bz_MTCoherogram(angle(hilbert(FiltFiltM(b,a,double(ca1.data(start:stop,1))))),angle(hilbert(FiltFiltM(b,a,double(ca3.data(start:stop,1))))),'range',[4 12],'window',.15,'step',.05);
            coh_ca1_ca3_phase(c2,:) = makeLength(mean(coherogram),1000);
            coh_ca1_ca3_z(c2,:) = zscore(coh_ca1_ca3(c2,:));
            coh_ca1_ca3_phase_z(c2,:) = zscore(coh_ca1_ca3_phase(c2,:));
            ca1_ca3_ind(c2) = ind;
            c2=c2+1;
            end
        end  
    end
    if ~isempty(ca3_ind)
    for i=1:3
        figure(i+10)
    subplot(3,2,1)
    plot(nanmean(coh_ca3_theta_gamma_z(ca3_ind==i,:)),'r')
    hold on
    plot(nanmean(coh_ca1_theta_gamma_z(ca1_ind==i,:)),'k')
    hold off
    subplot(3,2,2)
    plot(nanmean(coh_ca3_theta_gamma(ca3_ind==i,:)),'r')
    hold on
    plot(nanmean(coh_ca1_theta_gamma(ca1_ind==i,:)),'k')
    hold off
    subplot(3,2,3)
    imagesc(coh_ca3_theta_gamma_z(ca3_ind==i,:))
    subplot(3,2,4)
    imagesc(coh_ca1_theta_gamma_z(ca1_ind==i,:))
    subplot(3,2,5)
    plot(nanmean(coh_ca1_ca3_phase_z(ca1_ca3_ind==i,:)),'b')
    title('ca1/ca3 coherence')
    subplot(3,2,6)
    imagesc(coh_ca1_ca3_phase_z(ca1_ca3_ind==i,:))
    end
    pause(.01)
    end
%     end
%     end
    end
    clear ls ca1 ca3
    cd /home/david/datasets/lsDataset
%     cd /home/david/Dropbox/datasets/lsDataset
end