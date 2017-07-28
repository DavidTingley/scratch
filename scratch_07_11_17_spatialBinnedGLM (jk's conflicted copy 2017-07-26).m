% 
xml = LoadParameters;
warning off
load([xml.FileName '.behavior.mat'])
load([xml.FileName '.sessionInfo.mat'])
spikes = bz_GetSpikes;
lfp = bz_GetLFP(sessionInfo.thetaChans(2));
nBins = max(behavior.events.trials{1}.mapping);
% 
positionDecodingGLM_binnedspace_box.region = spikes.region;
positionDecodingGLM_binnedspace_box.sessionName = spikes.sessionName;
positionDecodingGLM_binnedspace_box.UID = spikes.UID;
for cell=1:length(spikes.times)
    positionDecodingGLM_binnedspace_box.results{cell} = table;
end
% % %    % set up phase coding data
[firingMaps] = bz_firingMap1D(spikes,behavior,lfp,4);
rateMap =firingMaps.rateMaps;
countMap =firingMaps.countMaps;
occuMap =firingMaps.occupancy;
phaseMap=firingMaps.phaseMaps;
[binnedPhaseMap] = bz_phaseMap2Bins(phaseMap,rateMap,behavior);
    
for smoothing = 1:round(nBins/2)
    disp(['smoothing by: ' num2str(smoothing) ' bins']);
    for cond = 1:length(unique(behavior.events.trialConditions))
%         figure(cond)
        % smooth data..
        for cell = 1:length(spikes.times)
                       %smoothing
            phase_trains_smooth=[];
            cos_phase_trains_smooth=[];
            sin_phase_trains_smooth=[];
            rates_trains_smooth = [];
            pos=[];
            for t = 1:size(binnedPhaseMap{cond},2)
                binnedPhaseMap{cond}(cell,t,isnan(binnedPhaseMap{cond}(cell,t,:)))=0;
                phase_trains_smooth=[phase_trains_smooth;...
                    circ_smoothTS(squeeze(binnedPhaseMap{cond}(cell,t,:)),smoothing,'method','mean','exclude',0)];
                
%                 phase_trains_smooth=[phase_trains_smooth;...
%                      Smooth(squeeze(binnedPhaseMap{cond}(cell,t,:)),smoothing,'type','c')];
%                 cos_phase_trains_smooth=[cos_phase_trains_smooth;...
%                      Smooth(squeeze(cos(binnedPhaseMap{cond}(cell,t,:))),smoothing,'type','c')];
%                 sin_phase_trains_smooth=[sin_phase_trains_smooth;...
%                      Smooth(squeeze(sin(binnedPhaseMap{cond}(cell,t,:))),smoothing,'type','c')];
                rates_trains_smooth = [rates_trains_smooth; ...
                                       smoothts(squeeze(countMap{cond}(cell,t,:))','b',smoothing)'];
                
                pos = [pos,[1:nBins]];                   
            end
%             subplot(2,1,1)
%             plot(rates_trains_smooth(1:200))
%             subplot(2,1,2)
%             plot(phase_trains_smooth(1:200))
%             title(smoothing)
%             pause

            r = rates_trains_smooth;
            p = phase_trains_smooth;
            p_cos =cos(phase_trains_smooth);
            p_sin = sin(phase_trains_smooth);
            
            count = 1;
            for iter = 1:10
            rr = randperm(length(r));
            pct = round(prctile(1:length(r),60));
            
            r_train = r(rr(1:pct));
            r_test = r(rr(pct+1:end));
            p_train = p(rr(1:pct));
            p_test = p(rr(pct+1:end));
            p_cos_train = p_cos(rr(1:pct));
            p_cos_test = p_cos(rr(pct+1:end));
            p_sin_train = p_sin(rr(1:pct));
            p_sin_test = p_sin(rr(pct+1:end));
            
            
            pos_train = pos(rr(1:pct));
            pos_test = pos(rr(pct+1:end));
            % rate
            [b dev stats] = glmfit([r_train],pos_train','normal');
            yfit = glmval(b,r_test,'identity');
            struct.mse_rate = nanmean((pos_test'-yfit).^2);
            struct.mse_rate_pval = stats.p';
            % phase all
            [b dev stats] = glmfit([p_train p_cos_train p_sin_train],pos_train','normal');
            yfit = glmval(b,[p_test p_cos_test p_sin_test],'identity');
            struct.mse_phase_all = nanmean((pos_test'-yfit).^2);
            struct.mse_phase_all_pval = stats.p';
            %phase
            [b dev stats] = glmfit([p_train],pos_train','normal');
            yfit = glmval(b,[p_test ],'identity');
            struct.mse_phase = nanmean((pos_test'-yfit).^2);
            struct.mse_phase_pval = stats.p';
            %phase cos
            [b dev stats] = glmfit([p_cos_train ],pos_train','normal');
            yfit = glmval(b,[ p_cos_test ],'identity');
            struct.mse_phase_cos = nanmean((pos_test'-yfit).^2);
            struct.mse_phase_cos_pval = stats.p';
            %phase sin
            [b dev stats] = glmfit([ p_sin_train],pos_train','normal');
            yfit = glmval(b,[ p_sin_test],'identity');
            struct.mse_phase_sin = nanmean((pos_test'-yfit).^2);
            struct.mse_phase_sin_pval = stats.p';
            % all predictors
            [b dev stats] = glmfit([r_train p_train p_cos_train p_sin_train],pos_train','normal');
            yfit = glmval(b,[r_test p_test p_cos_test p_sin_test],'identity');
            struct.mse_both  = nanmean((pos_test'-yfit).^2);
            struct.mse_both_pval = stats.p';
            
            [b dev stats] = glmfit([rand(length(pos_train),1)],pos_train','normal');
            yfit = glmval(b,[rand(length(pos_test),1)],'identity');
            struct.mse_chance  = nanmean((pos_test'-yfit).^2);
            
            % store peripherals
            struct.iter = count;
            struct.dfe = stats.dfe;
            struct.tau = smoothing;
            struct.condition = cond;
            positionDecodingGLM_binnedspace_box.results{cell} = [positionDecodingGLM_binnedspace_box.results{cell}; struct2table(struct)];
            
            count = 1+count;
            end
            if cell == 24 & cond == 1
%                 
                t_rate = varfun(@mean,positionDecodingGLM_binnedspace_box.results{cell},'InputVariables','mse_rate',...
                'GroupingVariables',{'tau','condition'});
                t_phase = varfun(@mean,positionDecodingGLM_binnedspace_box.results{cell},'InputVariables','mse_phase_all',...
                'GroupingVariables',{'tau','condition'});
                t_both = varfun(@mean,positionDecodingGLM_binnedspace_box.results{cell},'InputVariables','mse_both',...
                'GroupingVariables',{'tau','condition'});
                t_chance = varfun(@mean,positionDecodingGLM_binnedspace_box.results{cell},'InputVariables','mse_chance',...
                'GroupingVariables',{'tau','condition'});   
                t_std = varfun(@std,positionDecodingGLM_binnedspace_box.results{cell},'InputVariables','mse_chance',...
                'GroupingVariables',{'tau','condition'});  
                tab = join(join(join(t_rate,t_phase),t_both),t_chance);
                rows = find(tab.condition==cond);
                subplot(2,2,1);
%                 imagesc(phase_trains_smooth);
%                 caxis([-pi pi])
                subplot(2,2,2);
                boundedline(tab.tau(rows),tab.mean_mse_chance(rows),3*t_std.std_mse_chance(rows),'k')
                hold on
                plot(tab.tau(rows),tab.mean_mse_rate(rows),'r')
                plot(tab.tau(rows),tab.mean_mse_both(rows),'k')
                plot(tab.tau(rows),tab.mean_mse_phase_all(rows),'g')
                

%                 plot(positionDecodingGLM_binnedspace_box.results{cell}.tau(rows),...
%                     mean(positionDecodingGLM_binnedspace_box.results{cell}.mse_rate(rows,:)'),'r')
%                 hold on
%                 plot(positionDecodingGLM_binnedspace_box.results{cell}.tau(rows),...
%                     mean(positionDecodingGLM_binnedspace_box.results{cell}.mse_phase_all(rows,:)'),'g')
                plot(positionDecodingGLM_binnedspace_box.results{cell}.tau(rows),...
                    mean(positionDecodingGLM_binnedspace_box.results{cell}.mse_both(rows,:)'),'k')
                hold off
%                 subplot(2,2,3)
%                 imagesc(rates_trains_smooth);
%                 subplot(2,2,4)
                
                title([cell cond smoothing])
                pause(.1)
            end         
        end
        disp(['done with condition: ' num2str(cond) ' of ' num2str(length(unique(behavior.events.trialConditions)))]);
    end
    positionDecodingGLM_binnedspace_box.dateRun = date;  % this can take a very long time so lets save each loop...
% save([xml.FileName '.positionDecodingGLM_binnedspace_box.cellinfo.mat'],'positionDecodingGLM_binnedspace_box')
end
