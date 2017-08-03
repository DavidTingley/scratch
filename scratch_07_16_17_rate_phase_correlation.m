d  = dir('*201*');

hpc_corrs =[];
ls_corrs=[];
hpc_corrs_shuffle =[];
ls_corrs_shuffle=[];
hpc_corrs_change =[];
ls_corrs_change=[];
hpc_corrs_shuffle_change =[];
ls_corrs_shuffle_change=[];
ls_pvals = [];
hpc_pvals =[];
 
%% compile data

for ii=1:length(d)
   cd(d(ii).name)
   load([d(ii).name '.behavior.mat'])
    load([d(ii).name '.sessionInfo.mat'])
    spikes = bz_GetSpikes;
    load([d(ii).name '.firingMaps.cellinfo.mat'])
    
    rate_phase_corrs.UID = spikes.UID;
    rate_phase_corrs.sessionName = spikes.sessionName;
    rate_phase_corrs.region = spikes.region;
    
   for i=1:length(firingMaps.phaseMaps)
    for j=1:length(spikes.times)
        rate_phase_corrs.instRate_phase_corr{j,i}=NaN;
        rate_phase_corrs.instRate_phase_pval{j,i}=NaN;
        rate_phase_corrs.instRateChange_phaseChange_corr{j,i}=NaN;
        rate_phase_corrs.instRateChange_phaseChange_pval{j,i}=NaN;
        rate_phase_corrs.meanRateChangePhaseChangeCorr(j,i)=NaN;
        rate_phase_corrs.meanRatePhaseChangeCorr(j,i)=NaN;
        rate_phase_corrs.instRateChange_phaseChange_corr_shuffle{j,i}=NaN;
        rate_phase_corrs.instRateChange_phaseChange_pval_shuffle{j,i}=NaN;
        rate_phase_corrs.meanRateChangePhaseChangeCorr_shuffle(j,i)=NaN;
        rate_phase_corrs.instRate_phaseChange_corr_shuffle{j,i} = NaN; 
        rate_phase_corrs.instRate_phaseChange_pval_shuffle{j,i} = NaN; 
        rate_phase_corrs.meanRatePhaseChangeCorr_shuffle(j,i) = NaN; 
        rate_phase_corrs.ksPval(j,i) = NaN;
        
        if ~isempty(firingMaps.phaseMaps{i})
            if size(firingMaps.phaseMaps{i}{j},1) > 50  % min 5 spikes across all trials..
                for t = 1:size(firingMaps.rateMaps{i},2)
                    f = find(firingMaps.phaseMaps{i}{j}(:,2)==t);
                    if length(f) > 4  % min 3 spikes per trial...
                         [rate_phase_corrs.instRate_phase_corr{j,i}(t),rate_phase_corrs.instRate_phase_pval{j,i}(t)]= circ_corrcl(...
                            [firingMaps.phaseMaps{i}{j}(f,end)],firingMaps.phaseMaps{i}{j}(f,end-1));
                         [rate_phase_corrs.instRateChange_phaseChange_corr{j,i}(t),rate_phase_corrs.instRateChange_phaseChange_pval{j,i}(t)]= circ_corrcl(...
                            [firingMaps.phaseMaps{i}{j}(f,end)],[0; diff(firingMaps.phaseMaps{i}{j}(f,end-1))]);
                        for iter= 1:100
                            [rate_phase_corrs.instRate_phaseChange_corr_shuffle{j,i}(t,iter),rate_phase_corrs.instRate_phaseChange_pval_shuffle{j,i}(t,iter)]= circ_corrcl(...
                                [firingMaps.phaseMaps{i}{j}(f,end)],circshift([(firingMaps.phaseMaps{i}{j}(f,end-1))]',round(rand*100)));
                            [rate_phase_corrs.instRateChange_phaseChange_corr_shuffle{j,i}(t,iter),rate_phase_corrs.instRateChange_phaseChange_pval_shuffle{j,i}(t,iter)]= circ_corrcl(...
                                [firingMaps.phaseMaps{i}{j}(f,end)],circshift([0; diff(firingMaps.phaseMaps{i}{j}(f,end-1))]',round(rand*100)));
                        end
%                         if j == 80
%                             subplot(2,2,1)
%                             scatter([firingMaps.phaseMaps{i}{j}(f,end)],firingMaps.phaseMaps{i}{j}(f,end-1),'.k')
%                             hold on
%                             title('inst rate vs phase')
%                             subplot(2,2,2)
%                             scatter([firingMaps.phaseMaps{i}{j}(f,end)],[0; diff(firingMaps.phaseMaps{i}{j}(f,end-1))],'.k')
%                             hold on
%                             title('deriv inst rate vs phase')
%                             subplot(2,2,3)
%                             scatter((firingMaps.phaseMaps{i}{j}(f,end)),(firingMaps.phaseMaps{i}{j}(f,end-1)),'.k')
%                             hold on
%                             ylabel('delta rate')
%                             xlabel('delta phase')
%                         end
                   else  
                        for iter= 1:100
                            rate_phase_corrs.instRate_phaseChange_corr_shuffle{j,i}(t,iter)=nan;
                            rate_phase_corrs.instRate_phaseChange_pval_shuffle{j,i}(t,iter)= nan;
                            rate_phase_corrs.instRateChange_phaseChange_corr_shuffle{j,i}(t,iter)=nan;
                            rate_phase_corrs.instRateChange_phaseChange_pval_shuffle{j,i}(t,iter)= nan;
                        end
                         rate_phase_corrs.instRate_phase_corr{j,i}(t)=NaN;
                         rate_phase_corrs.instRate_phase_pval{j,i}(t)=NaN;
                         rate_phase_corrs.instRateChange_phaseChange_corr{j,i}(t)=NaN;
                         rate_phase_corrs.instRateChange_phaseChange_pval{j,i}(t)=NaN;
                    end
                     rate_phase_corrs.meanRateChangePhaseChangeCorr(j,i) = nanmean(rate_phase_corrs.instRateChange_phaseChange_corr{j,i});
                     rate_phase_corrs.meanRatePhaseChangeCorr(j,i) = nanmean(rate_phase_corrs.instRate_phase_corr{j,i});
                     rate_phase_corrs.meanRateChangePhaseChangeCorr_shuffle(j,i) = nanmean(rate_phase_corrs.instRateChange_phaseChange_corr_shuffle{j,i}(:));
                     rate_phase_corrs.meanRatePhaseChangeCorr_shuffle(j,i) = nanmean(rate_phase_corrs.instRate_phaseChange_corr_shuffle{j,i}(:));
                     if ~isnan(nanmean(rate_phase_corrs.instRate_phase_corr{j,i}))
                         [rate_phase_corrs.ksSig rate_phase_corrs.ksPval(j,i)] = ttest2(nanmean(rate_phase_corrs.instRate_phase_corr{j,i}),rate_phase_corrs.instRate_phaseChange_corr_shuffle{j,i}(:));
                     else
                         rate_phase_corrs.ksSig =nan;
                         rate_phase_corrs.ksPval(j,i)=nan;
                     end
                end
            end
        end
    end
   end
   for j=1:length(spikes.times)
       if strcmp(spikes.region{j},'hpc')
        subplot(4,2,1)
        hpc_corrs = [hpc_corrs ,rate_phase_corrs.meanRatePhaseChangeCorr(j,:)];
        hpc_corrs_shuffle = [hpc_corrs_shuffle ,rate_phase_corrs.meanRatePhaseChangeCorr_shuffle(j,:)];
        hpc_corrs_change = [hpc_corrs_change ,rate_phase_corrs.meanRateChangePhaseChangeCorr(j,:)];
        hpc_corrs_shuffle_change = [hpc_corrs_shuffle_change ,rate_phase_corrs.meanRateChangePhaseChangeCorr_shuffle(j,:)];
        hpc_pvals = [hpc_pvals,rate_phase_corrs.ksPval(j,:)];
        histogram(hpc_corrs,0:.02:1,'Normalization','pdf')
        hold on
        histogram(hpc_corrs_shuffle,0:.02:1,'Normalization','pdf','FaceColor','r')
        line([nanmedian(hpc_corrs(:)) nanmedian(hpc_corrs(:))],[0 20],'color','k')
        line([nanmedian(hpc_corrs_shuffle(:)) nanmedian(hpc_corrs_shuffle(:))],[0 20],'color','r')
        hold off
        %    set(gca,'yscale','log')
        title('hpc, phase vs inst rate corrs')
   
        subplot(4,2,2)
        histogram(hpc_corrs_change,0:.02:1,'Normalization','pdf')
        hold on
        histogram(hpc_corrs_shuffle_change,0:.02:1,'Normalization','pdf','FaceColor','r')
        line([nanmedian(hpc_corrs_change(:)) nanmedian(hpc_corrs_change(:))],[0 20],'color','k')
        line([nanmedian(hpc_corrs_shuffle_change(:)) nanmedian(hpc_corrs_shuffle_change(:))],[0 20],'color','r')
        hold off
        %    set(gca,'yscale','log')
        title('hpc, phase vs inst rate change corrs')
        
        subplot(4,2,3)
        scatter(hpc_corrs_change,hpc_corrs_shuffle_change,'.')
        subplot(4,2,4)
        histogram(hpc_pvals)
        
       elseif strcmp(spikes.region{j},'ls')
        subplot(4,2,5)
        ls_corrs = [ls_corrs ,rate_phase_corrs.meanRatePhaseChangeCorr(j,:)];
        ls_corrs_shuffle = [ls_corrs_shuffle ,rate_phase_corrs.meanRatePhaseChangeCorr_shuffle(j,:)];
        ls_corrs_change = [ls_corrs_change ,rate_phase_corrs.meanRateChangePhaseChangeCorr(j,:)];
        ls_corrs_shuffle_change = [ls_corrs_shuffle_change ,rate_phase_corrs.meanRateChangePhaseChangeCorr_shuffle(j,:)];
        ls_pvals = [ls_pvals,rate_phase_corrs.ksPval(j,:)];
        histogram(ls_corrs,0:.02:1,'Normalization','pdf')
        hold on
        histogram(ls_corrs_shuffle,0:.02:1,'Normalization','pdf','FaceColor','r')
        line([nanmedian(ls_corrs(:)) nanmedian(ls_corrs(:))],[0 20],'color','k')
        line([nanmedian(ls_corrs_shuffle(:)) nanmedian(ls_corrs_shuffle(:))],[0 20],'color','r')
        hold off
        title('ls, phase vs inst rate corrs')
        
        
        subplot(4,2,6)
        histogram(ls_corrs_change,0:.02:1,'Normalization','pdf')
        hold on
        histogram(ls_corrs_shuffle_change,0:.02:1,'Normalization','pdf','FaceColor','r')
        line([nanmedian(ls_corrs_change(:)) nanmedian(ls_corrs_change(:))],[0 20],'color','k')
        line([nanmedian(ls_corrs_shuffle_change(:)) nanmedian(ls_corrs_shuffle_change(:))],[0 20],'color','r')
        hold off
        title('ls, phase vs inst rate change corrs')
        
        subplot(4,2,7)
        scatter(ls_corrs_change,ls_corrs_shuffle_change,'.')
        subplot(4,2,8)
        histogram(ls_pvals)
       end
   end
   pause(.1)
   
   
%    save([rate_phase_corrs.sessionName '.rate_phase_corrs.cellinfo.mat'],'rate_phase_corrs')
   
clear rate_phase*
cd /home/david/datasets/lsDataset/
   
end