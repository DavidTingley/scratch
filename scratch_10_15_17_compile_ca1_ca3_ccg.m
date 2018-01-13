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
coh_ca3_ls_spk_pop = zeros(1,1000);
coh_ca1_ls_spk_pop = zeros(1,1000);
ca1_ca3_ind = [];
ccg_ca1_ca3 = zeros(1,1000,501);
ccg_ls_ca3 = zeros(1,1000,501);
ccg_ls_ca1 = zeros(1,1000,501);
ca1_recording=[];ca3_recording=[];
c1=1;
c2=1;
c3=1;
[b a] = butter(3,[4/625 12/625],'bandpass');

[bbb aaa] = butter(3,[50/625 500/625],'bandpass');
for ii=1:length(d)
    cd(d(ii).name)
    sessionInfo = bz_getSessionInfo;
    spikes = bz_GetSpikes('noprompts',true);
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
    end
    if ~isempty(sessionInfo.ls)
        ls = bz_GetLFP(sessionInfo.ls);
%         ls.data = FiltFiltM(b,a,double(ls.data));
    end
%     load([sessionInfo.FileName '.lfp.mat'],'ls','ca1','ca3')
% if ~isempty(sessionInfo.ca3)
%         for i=1:length(behavior.events.trialConditions)
%             if strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'linear')
%                 ind=1;
%             elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'central')
%                 ind=2;
%             elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'wheel')
%                 ind=3;
%             end
%             start = round(behavior.events.trialIntervals((i),1)*1250);
%             stop = round(behavior.events.trialIntervals((i),2)*1250);
%             ccc = zeros(1000,501);
%             if stop-start > 1000
%             for step = 1:stop-start
%                     [cc(step,:)] = crosscorr(cos(angle(hilbert(FiltFiltM(b,a,double(ls.data(start+step-250:start+step+250,1)))))),cos(angle(hilbert(FiltFiltM(b,a,double(ca3.data(start+step-250:start+step+250,1)))))),250);
%             end  
%             for step=1:501
%                 ccc(:,step) = makeLength(cc(:,step),1000);
%             end
%             end
%             ccg_ls_ca1(c1,:,:) = ccc; clear cc ccc
%             c1=c1+1;
%             end
%         end  
%     if ~isempty(sessionInfo.ca1)
%         for i=1:length(behavior.events.trialConditions)
%             if strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'linear')
%                 ind=1;
%             elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'central')
%                 ind=2;
%             elseif strcmp(behavior.events.conditionType{behavior.events.trialConditions(i)},'wheel')
%                 ind=3;
%             end
%             start = round(behavior.events.trialIntervals((i),1)*1250);
%             stop = round(behavior.events.trialIntervals((i),2)*1250);
%             ccc = zeros(1000,501);
%             if stop-start > 1000
%             for step = 1:stop-start
%                     [cc(step,:)] = crosscorr(cos(angle(hilbert(FiltFiltM(b,a,double(ls.data(start+step-250:start+step+250,1)))))),cos(angle(hilbert(FiltFiltM(b,a,double(ca1.data(start+step-250:start+step+250,1)))))),250);
%             end  
%             for step=1:501
%                 ccc(:,step) = makeLength(cc(:,step),1000);
%             end
%             end
%             ccg_ls_ca3(c2,:,:) = ccc; clear cc ccc
%             c2=c2+1;
%             end
%         end  
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
            ccc = zeros(1000,501);
            if stop-start > 1000
                ccg_ca1_ca3(c3,:) = makeLength(circ_dist(angle(hilbert(FiltFiltM(b,a,double(ca1.data(start:stop,1))))),angle(hilbert(FiltFiltM(b,a,double(ca3.data(start:stop,1)))))),1000);
%             for step = 1:stop-start
%                     [cc(step,:)] = crosscorr(cos(angle(hilbert(FiltFiltM(b,a,double(ca1.data(start+step-250:start+step+250,1)))))),cos(angle(hilbert(FiltFiltM(b,a,double(ca3.data(start+step-250:start+step+250,1)))))),250);
%             end  
%             for step=1:501
%                 ccc(:,step) = makeLength(cc(:,step),1000);
%             end
            end
%             ccg_ca1_ca3(c3,:,:) = ccc; clear cc ccc
            c3=c3+1;
            end
        end  
    end
    
    for i=1:3
        figure(i+10)
    subplot(3,2,1)
    imagesc(squeeze(nanmean(ccg_ca1_ca3))')
    subplot(3,2,2)
    plot(nanmean(squeeze(nanmean(ccg_ca1_ca3))'))
    subplot(3,2,3)
    imagesc(squeeze(nanmean(ccg_ls_ca3))')
    subplot(3,2,4)
    plot(nanmean(squeeze(nanmean(ccg_ls_ca3))'))
    subplot(3,2,5)
    imagesc(squeeze(nanmean(ccg_ls_ca1))')
    subplot(3,2,6) 
    plot(nanmean(squeeze(nanmean(ccg_ls_ca1))'))
    pause(.01)
%     end
%     end
    end
    clear ls ca1 ca3 spk 
    cd /home/david/datasets/lsDataset
%     cd /home/david/Dropbox/datasets/lsDataset
end