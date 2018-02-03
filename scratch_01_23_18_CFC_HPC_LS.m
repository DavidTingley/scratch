% clear 
d = dir('*201*');
% clf
[b a] = butter(4,[6/625 10/625],'bandpass');

 ord = randperm(length(d));

for rec= 51:length(d)
    rec
    cd(d((rec)).name)
    sessionInfo = bz_getSessionInfo;

    if ~isempty(sessionInfo.ca1)
        hpc = bz_GetLFP(sessionInfo.ca1);
        label(rec) = 1;
    else
        hpc = bz_GetLFP(sessionInfo.ca3);
        label(rec) = 3;
    end

    ls = bz_GetLFP(sessionInfo.ls);
    load([sessionInfo.FileName '.behavior.mat'],'behavior')
    
    for trial = 1:length(behavior.events.trialIntervals)
    [blah start] = min(abs(hpc.timestamps-min(behavior.events.trialIntervals(trial,1))));
    [blah stop] = min(abs(hpc.timestamps-max(behavior.events.trialIntervals(trial,2))));
    freqs2=25:2:200;
    freqs1=5:0.5:10;
    if stop-start>2500
    [cmats_sig_hpc_ls{rec}(trial,:,:,:,:),strength,pref_phase,mean_phase] = thetamod2(double(hpc.data(start:stop)),double(ls.data(start:stop)),freqs2,[2 3 5 6],1250,[25 200]);
    [tort_hpc_ls{rec}(trial,:,:),freq,MI] = CFCtort(double(hpc.data(start:stop)),double(ls.data(start:stop)),freqs1,freqs2,1250,[30 200]);
    
    [cmats_sig_hpc_hpc{rec}(trial,:,:,:,:),strength,pref_phase,mean_phase] = thetamod2(double(hpc.data(start:stop)),double(hpc.data(start:stop)),freqs2,[2 3 5 6],1250,[25 200]);
    [tort_hpc_hpc{rec}(trial,:,:),freq,MI] = CFCtort(double(hpc.data(start:stop)),double(hpc.data(start:stop)),freqs1,freqs2,1250,[30 200]);
    
    [cmats_sig_ls_ls{rec}(trial,:,:,:,:),strength,pref_phase,mean_phase] = thetamod2(double(ls.data(start:stop)),double(ls.data(start:stop)),freqs2,[2 3 5 6],1250,[25 200]);
    [tort_ls_ls{rec}(trial,:,:),freq,MI] = CFCtort(double(ls.data(start:stop)),double(ls.data(start:stop)),freqs1,freqs2,1250,[30 200]);
    end
    end
%     
%     %% HPC/LS coupling
%     for j =1:629
%        f = find(hpc_phases==j-315); 
%        hpc_ls_cfc((rec),j,:) = nanmean(ls_power(:,f),2)'; 
%     end
%     for j=1:50
%         hpc_ls_cfc((rec),:,j) = zscore(hpc_ls_cfc((rec),:,j));
%     end
% %     hpc_ls_powpow((rec),:,:) = corr(hpc_power',ls_power');
% 
%     %% HPC/HPC coupling
%     for j =1:629
%        f = find(hpc_phases==j-315); 
%        hpc_hpc_cfc((rec),j,:) = nanmean(hpc_power(:,f),2)'; 
%     end
%     for j=1:50
%         hpc_hpc_cfc((rec),:,j) = zscore(hpc_hpc_cfc((rec),:,j));
%     end
% %     hpc_hpc_powpow((rec),:,:) = corr(hpc_power',hpc_power');
% 
%     %% LS/LS coupling
%     for j =1:629
%        f = find(ls_phases==j-315); 
%        ls_ls_cfc((rec),j,:) = nanmean(ls_power(:,f),2)'; 
%     end
%     for j=1:50
%         ls_ls_cfc((rec),:,j) = zscore(ls_ls_cfc((rec),:,j));
%     end
%     ls_ls_powpow((rec),:,:) = corr(ls_power',ls_power');
            
    %% plotting
%     subplot(3,2,1)
%     imagesc(-3.14:6.29/100:3.14,freqs,squeeze(nanmean(hpc_ls_cfc(label==1,:,:),1))');
%     title('CA1/ls')
%     
%     subplot(3,2,2)
%     imagesc(-3.14:6.29/100:3.14,freqs,squeeze(nanmean(hpc_ls_cfc(label==3,:,:),1))');
%     title('CA3/ls')
% %     subplot(3,2,2)
% %     imagesc(freqs,freqs,hpc_ls_powpow)
%     
%     subplot(3,2,3)
%     imagesc(-3.14:6.29/100:3.14,freqs,squeeze(nanmean(hpc_hpc_cfc,1))');
%     title('hpc/hpc')
% %     subplot(3,2,4)
% %     imagesc(freqs,freqs,hpc_hpc_powpow)
%     
%     subplot(3,2,5)
%     imagesc(-3.14:6.29/100:3.14,freqs,squeeze(nanmean(ls_ls_cfc,1))');
%     title('ls/ls')
%     subplot(3,2,6)
%     imagesc(freqs,freqs,ls_ls_powpow)
    
    pause(.01)
%         subplot(3,2,1)
%     imagesc(squeeze(mean(hpc_ls_cfc)));
    cd /home/david/datasets/lsDataset/
    rec
%     save('/home/david/Dropbox/HPC_LS_CFC_dataset.mat','-v7.3')
end

