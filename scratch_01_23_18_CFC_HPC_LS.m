% clear 
d = dir('*201*');
% clf
[b a] = butter(3,[6/625 10/625],'bandpass');

ord = randperm(length(d));
cfc_hpc_hpc = zeros(68,170,1001);
cfc_hpc_ls = zeros(68,170,1001);
cfc_ls_ls = zeros(68,170,1001);
    
for rec= 1:length(d)
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
    
    [pks locs]= findpeaks(-FiltFiltM(b,a,double(hpc.data(start:stop))));
    [freqs,t,hpc_spec] = bz_WaveSpec(double(hpc.data(start:stop)),[30 200],170,3,1/1250,'lin');
    hpc_spec = abs(hpc_spec);
    
    [freqs,t,ls_spec] = bz_WaveSpec(double(ls.data(start:stop)),[30 200],170,3,1/1250,'lin');
    ls_spec = abs(ls_spec);
    
    for pk = 1:length(pks)
        if locs(pk) > 500 & locs(pk) < stop-start - 500
            cfc_hpc_hpc(rec,:,:) = squeeze(cfc_hpc_hpc(rec,:,:)) + hpc_spec(:,locs(pk)-500:locs(pk)+500);
            cfc_hpc_ls(rec,:,:) = squeeze(cfc_hpc_ls(rec,:,:)) + ls_spec(:,locs(pk)-500:locs(pk)+500);
        end
    end
    [pks locs]= findpeaks(-FiltFiltM(b,a,double(ls.data(start:stop))));
    for pk = 1:length(pks)
        if locs(pk) > 500 & locs(pk) < stop-start -500
            cfc_ls_ls(rec,:,:) = squeeze(cfc_ls_ls(rec,:,:)) + ls_spec(:,locs(pk)-500:locs(pk)+500);
        end
    end
%     freqs2=25:2:400;
%     freqs1=5:0.5:10;
%     if stop-start>5000
%     [cmats_sig_hpc_ls{rec}(trial,:,:,:,:),strength,pref_phase,mean_phase] = thetamod2(double(hpc.data(start:stop)),double(ls.data(start:stop)),freqs2,[5 10],1500,[25 400]);
%     [tort_hpc_ls{rec}(trial,:,:),freq,MI] = CFCtort(double(hpc.data(start:stop)),double(ls.data(start:stop)),freqs1,freqs2,1500,[30 400]);
%     
%     [cmats_sig_hpc_hpc{rec}(trial,:,:,:,:),strength,pref_phase,mean_phase] = thetamod2(double(hpc.data(start:stop)),double(hpc.data(start:stop)),freqs2,[5 10],1500,[25 400]);
%     [tort_hpc_hpc{rec}(trial,:,:),freq,MI] = CFCtort(double(hpc.data(start:stop)),double(hpc.data(start:stop)),freqs1,freqs2,1500,[30 400]);
%     
%     [cmats_sig_ls_ls{rec}(trial,:,:,:,:),strength,pref_phase,mean_phase] = thetamod2(double(ls.data(start:stop)),double(ls.data(start:stop)),freqs2,[5 10],1500,[25 400]);
%     [tort_ls_ls{rec}(trial,:,:),freq,MI] = CFCtort(double(ls.data(start:stop)),double(ls.data(start:stop)),freqs1,freqs2,1500,[30 400]);
%     
%     else
%         cmats_sig_hpc_ls{rec}(trial,:,:,:,:) = nan(88,32,1,2);
%         tort_hpc_ls{rec}(trial,:,:) = nan(11,88);
%         cmats_sig_ls_ls{rec}(trial,:,:,:,:) = nan(88,32,1,2);
%         tort_ls_ls{rec}(trial,:,:) = nan(11,88);
%         cmats_sig_hpc_hpc{rec}(trial,:,:,:,:) = nan(88,32,1,2);
%         tort_hpc_hpc{rec}(trial,:,:) = nan(11,88);
%         
%     end
    end

    %% plotting

%     nbins_theta = 32;
%     xtheta = 0:(2*pi/nbins_theta):(2*pi);
    subplot(3,2,1)
    imagesc(-400:400,30:200,squeeze(nanmean(cfc_hpc_hpc)))
    for i=1:170
        cfc_hpc_hpc_z(i,:) = zscore(nanmean(cfc_hpc_hpc(:,i,:)));
    end
    subplot(3,2,2)
    imagesc(-400:400,30:200,cfc_hpc_hpc_z)
%     mat_hpc_ls(rec,:,:) = squeeze(nanmean(cmats_sig_hpc_ls{rec}(:,:,:,1,1),1));
%     imagesc([midpoints(xtheta(1:2)) 4*pi-midpoints(xtheta(1:2))]*180/pi, freqs2, squeeze(nanmean(mat_hpc_ls(:,:,[1:end 1:end]),1)))
%     hold on
%     plot([midpoints(xtheta) midpoints(xtheta)+2*pi]*180/pi, (cos([midpoints(xtheta) midpoints(xtheta)]) - 1)*10 + 150, 'k--');
    title('hpc/hpc')
    colorbar
    
    subplot(3,2,3)
    imagesc(-400:400,30:200,squeeze(nanmean(cfc_hpc_ls)))
    for i=1:170
        cfc_hpc_ls_z(i,:) = zscore(nanmean(cfc_hpc_ls(:,i,:)));
    end
    subplot(3,2,4)
    imagesc(-400:400,30:200,cfc_hpc_ls_z)
%     mat_hpc_hpc(rec,:,:) = squeeze(nanmean(cmats_sig_hpc_hpc{rec}(:,:,:,1,1),1));
%     imagesc([midpoints(xtheta(1:2)) 4*pi-midpoints(xtheta(1:2))]*180/pi, freqs2, squeeze(nanmean(mat_hpc_hpc(:,:,[1:end 1:end]),1)))
%     hold on
%     plot([midpoints(xtheta) midpoints(xtheta)+2*pi]*180/pi, (cos([midpoints(xtheta) midpoints(xtheta)]) - 1)*10 + 150, 'k--');
    title('hpc/ls')
    colorbar
    
    subplot(3,2,5)
    imagesc(-400:400,30:200,squeeze(nanmean(cfc_ls_ls)))
    for i=1:170
        cfc_ls_ls_z(i,:) = zscore(nanmean(cfc_ls_ls(:,i,:)));
    end
    subplot(3,2,6)
    imagesc(-400:400,30:200,cfc_ls_ls_z)
%     mat_ls_ls(rec,:,:) = squeeze(nanmean(cmats_sig_ls_ls{rec}(:,:,:,1,1),1));
%     imagesc([midpoints(xtheta(1:2)) 4*pi-midpoints(xtheta(1:2))]*180/pi, freqs2, squeeze(nanmean(mat_ls_ls(:,:,[1:end 1:end]),1)))
%     hold on
%     plot([midpoints(xtheta) midpoints(xtheta)+2*pi]*180/pi, (cos([midpoints(xtheta) midpoints(xtheta)]) - 1)*10 + 150, 'k--');
    title('ls/ls')
    colorbar
    
    
    pause(.01)
    cd /home/david/datasets/lsDataset/
    rec
    save('/home/david/Dropbox/HPC_LS_CFC_dataset.mat','-v7.3')
end

