

% 
% xml = LoadParameters;
% 
% load([xml.FileName '.behavior.mat'])
% load([xml.FileName '.sessionInfo.mat'])
% spikes = bz_GetSpikes;
% lfp = bz_GetLFP(sessionInfo.thetaChans(2));
% nBins = max(behavior.events.trials{1}.mapping);
% 
positionDecodingGLM_binnedspace.region = spikes.region;
positionDecodingGLM_binnedspace.sessionName = spikes.sessionName;
positionDecodingGLM_binnedspace.UID = spikes.UID;
for cell=1:length(spikes.times)
    positionDecodingGLM_binnedspace.results{cell} = table;
end
% % %    % set up phase coding data
% [rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,4);
% [binnedPhaseMap] = bz_phaseMap2Bins(phaseMap,rateMap,behavior);
    
for smoothing = 1:round(nBins/2)
    disp(['smoothing by: ' num2str(smoothing) ' bins']);
    for cond = 5% 1:length(unique(behavior.events.trialConditions))
%         figure(cond)
        % smooth data..
        for cell = 1:length(spikes.times)
           for trial = 1:size(binnedPhaseMap{cond},2)
              binnedPhaseMap_smooth{cond}(cell,trial,:) = circ_smoothTS(squeeze(binnedPhaseMap{cond}(cell,trial,:)),smoothing,'method','mean','exclude',0); 
              binnedPhaseMap_smooth{cond}(cell,trial,binnedPhaseMap_smooth{cond}(cell,trial,:)==0)=nan;
              rateMap_smooth{cond}(cell,trial,:) = smooth(squeeze(countMap{cond}(cell,trial,:))',smoothing);
           end
           range = 0:max(max(squeeze(rateMap_smooth{cond}(cell,:,:))))./63:max(max(squeeze(rateMap_smooth{cond}(cell,:,:))));
           if isempty(range)
               range = [0 1];
           end
        end
        binnedPhaseMap_smooth{cond}(isnan(binnedPhaseMap_smooth{cond}))=0;    
        f = find(rateMap_smooth{cond}(cell,:)==0);
        rateMap_smooth{cond}(cell,f)=nan;
        f = find(binnedPhaseMap_smooth{cond}(cell,:)==0);
        binnedPhaseMap_smooth{cond}(cell,f)=nan;
%         rateMap_disc{cond}(cell,:,:) = discretize(rateMap_smooth{cond}(cell,:,:),range);  % discretize both rate/phase to same # of bins...
%         phaseMap_disc{cond}(cell,:,:) = discretize(binnedPhaseMap_smooth{cond}(cell,:,:),-pi:.1:pi);% ,-1:.032:1);
%         phaseMap_disc_sin{cond}(cell,:,:) = discretize(sin(binnedPhaseMap_smooth{cond}(cell,:,:)),-pi:.1:pi);
%         phaseMap_disc_cos{cond}(cell,:,:) = discretize(cos(binnedPhaseMap_smooth{cond}(cell,:,:)),-pi:.1:pi);
        % compile data
        for cell = 80 %1:length(spikes.times)
            r = squeeze(rateMap_smooth{cond}(cell,:,:));
            p = squeeze(binnedPhaseMap_smooth{cond}(cell,:,:));
            z = repmat([1:nBins]',size(r,1),1);
            r = reshape(r',size(r,1)*size(r,2),1);
            p = reshape(p',size(p,1)*size(p,2),1);
            
            count = 1;
            for iter = 1:10
            rr = randperm(length(r));
            pct = round(prctile(1:length(r),60));
            r_train = r(rr(1:pct));
            r_test = r(rr(pct+1:end));
            p_train = p(rr(1:pct));
            p_test = p(rr(pct+1:end));
            z_train = z(rr(1:pct));
            z_test = z(rr(pct+1:end));
            
            [b dev stats] = glmfit([r_train],z_train,'normal');
            yfit = glmval(b,r_test,'identity');
            struct.mse_rate(count) = nanmean((z_test-yfit).^2);
%             struct.mse_rate_pval(count,:) = stats.p';
            
            [b dev stats] = glmfit([p_train cos(p_train) sin(p_train)],z_train,'normal');
            yfit = glmval(b,[p_test cos(p_test) sin(p_test)],'identity');
            struct.mse_phase_all(count) = nanmean((z_test-yfit).^2);
%             struct.mse_phase_all_pval(count,:) = stats.p';
            
            [b dev stats] = glmfit([p_train],z_train,'normal');
            yfit = glmval(b,[p_test ],'identity');
            struct.mse_phase(count) = nanmean((z_test-yfit).^2);
%             struct.mse_phase_pval(count,:) = stats.p';
            
            [b dev stats] = glmfit([cos(p_train) ],z_train,'normal');
            yfit = glmval(b,[ cos(p_test) ],'identity');
            struct.mse_phase_cos(count) = nanmean((z_test-yfit).^2);
%             struct.mse_phase_cos_pval(count,:) = stats.p';
            
            [b dev stats] = glmfit([ sin(p_train)],z_train,'normal');
            yfit = glmval(b,[ sin(p_test)],'identity');
            struct.mse_phase_sin(count) = nanmean((z_test-yfit).^2);
%             struct.mse_phase_sin_pval(count,:) = stats.p';
            
            [b dev stats] = glmfit([r_train p_train cos(p_train) sin(p_train)],z_train,'normal');
            yfit = glmval(b,[r_test p_test cos(p_test) sin(p_test)],'identity');
            struct.mse_both(count)  = nanmean((z_test-yfit).^2);
%             struct.mse_both_pval(count,:) = stats.p';

            % reshape data
%             struct.mse_both = struct.mse_both';
%             struct.mse_rate = struct.mse_rate';
%             struct.mse_phase = struct.mse_phase';
%             struct.mse_phase_cos = struct.mse_phase_cos';
%             struct.mse_phase_sin = struct.mse_phase_sin';
%             struct.mse_phase_all = struct.mse_phase_all';
            
            % store peripherals
            struct.iter = count;
            struct.dfe = stats.dfe;
            struct.tau = smoothing;
            struct.condition = cond;
            positionDecodingGLM_binnedspace.results{cell} = [positionDecodingGLM_binnedspace.results{cell}; struct2table(struct)];
            
            count = 1+count;
            end
            if cell == 80 & cond == 5
                rows = find(positionDecodingGLM_binnedspace.results{cell}.condition==cond);
                subplot(2,2,1);
                imagesc((squeeze(binnedPhaseMap_smooth{cond}(cell,:,:))));
                caxis([-1 1])
                subplot(2,2,2);
                plot(positionDecodingGLM_binnedspace.results{cell}.tau(rows),...
                    mean(positionDecodingGLM_binnedspace.results{cell}.mse_rate(rows,:)'),'r')
                hold on
                plot(positionDecodingGLM_binnedspace.results{cell}.tau(rows),...
                    mean(positionDecodingGLM_binnedspace.results{cell}.mse_phase_all(rows,:)'),'g')
                plot(positionDecodingGLM_binnedspace.results{cell}.tau(rows),...
                    mean(positionDecodingGLM_binnedspace.results{cell}.mse_both(rows,:)'),'k')
                hold off
                subplot(2,2,3)
                imagesc(squeeze(rateMap_smooth{cond}(cell,:,:)));
%                 subplot(2,2,4)
%                 scatter(phaseMap{cond}{cell}(:,1),phaseMap{cond}{cell}(:,end)+2*pi,'.k');
%                 title([cell cond smoothing])
                pause(.1)
            end         
        end
        disp(['done with condition: ' num2str(cond) ' of ' num2str(length(unique(behavior.events.trialConditions)))]);
    end
    positionDecodingGLM_binnedspace.dateRun = date;  % this can take a very long time so lets save each loop...
% save([xml.FileName '.positionDecodingGLM_binnedspace.cellinfo.mat'],'positionDecodingGLM_binnedspace')
end

