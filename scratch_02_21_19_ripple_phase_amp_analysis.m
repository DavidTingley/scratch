addpath(genpath('/home/david/Dropbox/code'))
cd /home/david/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
% clear all
% phaseAmp = zeros(60,49);
warning off

d = dir('*201*');

for rec = 1:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    sessionInfo.FileName;
    ls_ripples = bz_LoadEvents(pwd,'LSRipples');
    ca1_ripples = bz_LoadEvents(pwd,'CA1Ripples');
    sleep = bz_LoadStates(pwd,'SleepState');

    if ~isempty(ls_ripples) & ~isempty(sessionInfo.ls)
        lfp = bz_GetLFP(sessionInfo.ls);
        pow_corr_t = nan(length(ls_ripples.peaks),201,201);
        pow_corr_shuf_t = nan(length(ls_ripples.peaks),1,201,201);
        
%         parfor i=1:length(ls_ripples.peaks)
%             i
            l = lfp;
%             [a start] = min(abs(ls_ripples.timestamps(i,1)-l.timestamps));
%             [a stop] = min(abs(ls_ripples.timestamps(i,2)-l.timestamps));
%             if start > 2500 & stop < length(l.data)-2500
%             l.data = l.data(start-2500:stop+2500);
%             l.timestamps = l.timestamps(start-2500:stop+2500);
%             [pow_corr_t(i,:,:) pow_corr_shuf_t(i,:,:,:) phaseAmps(i,:,:)] = runCorr(l,ls_ripples);
            [pow_corr_t(i,:,:) pow_corr_shuf_t(i,:,:,:) ] = runCorr(l,ls_ripples);
%             end
%         end
        pow_corr{rec} = pow_corr_t; clear pow_corr_t;
        pow_corr_shuf{rec} = pow_corr_shuf_t; clear pow_corr_shuf_t;
%         phaseAmp = squeeze(nansum(phaseAmps)); clear phaseAmps;
        
    subplot(2,2,1)
    p_co(rec,:,:) = squeeze(nanmean(pow_corr{rec}));
    imagesc(squeeze(nanmean(p_co)))
    subplot(2,2,2)
%     imagesc([5:50],[120:250],squeeze((phaseAmp)))
    subplot(2,2,3)
    p_co_shuf(rec,:,:) = squeeze(nanmean(nanmean(pow_corr_shuf{rec},2)));
    imagesc(squeeze(nanmean(p_co_shuf)))
    
    pause(.1);
    end
    save('/home/david/Dropbox/LS_ripple_phase_amp_analysis_long.mat','-v7.3')
   cd /home/david/datasets/ripples_LS 
% cd E:\datasets\ripples_LS
   clear coupling
end

function [pow_corr pow_corr_shuf ] = runCorr(l,ls_ripples)

        pow = zeros(201,length(l.data),'single');
        parfor ii=1:201
%                 [b,a] = butter(4,[(ii)/625 (ii+3)/625],'bandpass');
            pow(ii,:) = single(fastrms(bz_Filter(double(l.data),'order',...
                3,'filter','butter','passband',[ii ii+1]),round(2000/ii)));
            ii
        end
        for ii=1:201
            for j=i:201
                pow_corr(ii,j) = corr(pow(ii,:)',pow(j,:)');
                pow_corr(j,ii) = pow_corr(ii,j); 
                for iter = 1
                    pow_corr_shuf(iter,ii,j) = corr(circshift(pow(ii,:)',randi([1 length(pow)],1)),pow(j,:)');
                end
            end
        end
        clear pow;

%         phaseAmp = bz_ModIndex(l,[5:50],[120:250],0);

end