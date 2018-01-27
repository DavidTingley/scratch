clear 
d = dir('*201*');
clf
[b a] = butter(4,[4/625 10/625],'bandpass');

ord = randperm(length(d));

for rec= 1:length(d)
    cd(d(ord(rec)).name)
    sessionInfo = bz_getSessionInfo;

    if ~isempty(sessionInfo.ca1)
        hpc = bz_GetLFP(sessionInfo.ca1);
    else
        hpc = bz_GetLFP(sessionInfo.ca3);
    end

    ls = bz_GetLFP(sessionInfo.ls);
    load([sessionInfo.FileName '.behavior.mat'],'behavior')
    
    
    [blah start] = min(abs(hpc.timestamps-behavior.timestamps(1)));
    [blah stop] = min(abs(hpc.timestamps-behavior.timestamps(end)));
    
    hpc_phases = angle(hilbert(FiltFiltM(b,a,double(hpc.data(start:stop)))));
    ls_phases = angle(hilbert(FiltFiltM(b,a,double(ls.data(start:stop)))));
    hpc_phases = round(hpc_phases*100);
    ls_phases = round(ls_phases*100);
    
    [freqs,t,spec] = WaveSpec(hpc.data(start:stop),[1 200],200,3,1/1250,'lin');
    hpc_power = single(abs(spec));
    [freqs,t,spec] = WaveSpec(ls.data(start:stop),[1 200],200,3,1/1250,'lin');
    ls_power = single(abs(spec));

    hpc_phases(mean(hpc_power(4:10,:))<mean(mean(hpc_power(4:10,:)))+3*std(mean(hpc_power(4:10,:)))) = nan;
    ls_phases(mean(ls_power(4:10,:))<mean(mean(ls_power(4:10,:)))+3*std(mean(ls_power(4:10,:)))) = nan;
    
    hpc_power(:,mean(hpc_power(4:10,:))<mean(mean(hpc_power(4:10,:)))+3*std(mean(hpc_power(4:10,:)))) = nan;
    ls_power(:,mean(ls_power(4:10,:))<mean(mean(ls_power(4:10,:)))+3*std(mean(ls_power(4:10,:)))) = nan;
    
    clear spec;
    
    %% HPC/LS coupling
    for j =1:629
       f = find(hpc_phases==j-315); 
       hpc_ls_cfc(ord(rec),j,:) = nanmean(ls_power(:,f),2)'; 
    end
    for j=1:200
        hpc_ls_cfc(ord(rec),:,j) = minmax_norm(hpc_ls_cfc(ord(rec),:,j));
    end
%     hpc_ls_powpow(ord(rec),:,:) = corr(hpc_power',ls_power');

    %% HPC/HPC coupling
    for j =1:629
       f = find(hpc_phases==j-315); 
       hpc_hpc_cfc(ord(rec),j,:) = nanmean(hpc_power(:,f),2)'; 
    end
    for j=1:200
        hpc_hpc_cfc(ord(rec),:,j) = minmax_norm(hpc_hpc_cfc(ord(rec),:,j));
    end
%     hpc_hpc_powpow(ord(rec),:,:) = corr(hpc_power',hpc_power');

    %% LS/LS coupling
    for j =1:629
       f = find(ls_phases==j-315); 
       ls_ls_cfc(ord(rec),j,:) = nanmean(ls_power(:,f),2)'; 
    end
    for j=1:200
        ls_ls_cfc(ord(rec),:,j) = minmax_norm(ls_ls_cfc(ord(rec),:,j));
    end
%     ls_ls_powpow(ord(rec),:,:) = corr(ls_power',ls_power');
            
    %% plotting
    subplot(3,2,1)
    imagesc(-3.14:6.29/100:3.14,freqs,squeeze(nanmean(hpc_ls_cfc,1))');
    title('hpc/ls')
%     subplot(3,2,2)
%     imagesc(freqs,freqs,hpc_ls_powpow)
    
    subplot(3,2,3)
    imagesc(-3.14:6.29/100:3.14,freqs,squeeze(nanmean(hpc_hpc_cfc,1))');
    title('hpc/hpc')
%     subplot(3,2,4)
%     imagesc(freqs,freqs,hpc_hpc_powpow)
    
    subplot(3,2,5)
    imagesc(-3.14:6.29/100:3.14,freqs,squeeze(nanmean(ls_ls_cfc,1))');
    title('ls/ls')
%     subplot(3,2,6)
%     imagesc(freqs,freqs,ls_ls_powpow)
    
    pause(.01)
%         subplot(3,2,1)
%     imagesc(squeeze(mean(hpc_ls_cfc)));
    cd /home/david/datasets/lsDataset/
    save('/home/david/Dropbox/HPC_LS_CFC_dataset.mat','-v7.3')
end

