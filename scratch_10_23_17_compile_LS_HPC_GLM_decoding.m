d = dir('*201*');
tau = 35;
count=1;

for ii=1:length(d)
    cd(d(ii).name)
    disp(['working on ' d(ii).name ''])
    spikes = bz_GetSpikes;
    sessionInfo = bz_getSessionInfo;
    load([sessionInfo.FileName '.behavior.mat'],'behavior')
    load([sessionInfo.FileName '.positionDecodingGLM_binnedspace_box.cellinfo.mat'])
    
    for cell=1:length(spikes.times)
        if strcmp(spikes.region{cell},'ls')
            if exist([sessionInfo.FileName '_rate_decoding_GLM_cell_' num2str(cell) '.mat'])
            load([sessionInfo.FileName '_rate_decoding_GLM_cell_' num2str(cell) '.mat'],'lsRateDecodingHPC_POP_MaxCorr')
            if ~isempty(lsRateDecodingHPC_POP_MaxCorr.results)
            conditions = unique(lsRateDecodingHPC_POP_MaxCorr.results.condition);
            for cond = 1:length(conditions)
                if isempty(strmatch('fits',lsRateDecodingHPC_POP_MaxCorr.results.Properties.VariableNames))
                    warning(['missing .fits for ' sessionInfo.FileName ', cell #' num2str(cell)])
                else
                    % check if LS cell precesses...
                    rows = find(positionDecodingGLM_binnedspace_box.results{cell}.condition==cond);
                    cols = find(positionDecodingGLM_binnedspace_box.results{cell}.tau==80);
                    rows = intersect(rows,cols);
                    if ~isempty(rows)
                    precess_phase = positionDecodingGLM_binnedspace_box.results{cell}.mse_phase_all(rows);
                    precess_chance = positionDecodingGLM_binnedspace_box.results{cell}.mse_chance(rows);
                     precess_rate = positionDecodingGLM_binnedspace_box.results{cell}.mse_rate(rows);
                    [a pval] = kstest2(precess_phase,precess_chance);
                    [a pval_r] = kstest2(precess_rate,precess_chance);
%                     if pval < .01
%                         if strcmp(behavior.events.conditionType{cond},'wheel') 
                            rows = find(lsRateDecodingHPC_POP_MaxCorr.results.tau==tau);
                            cols = find(lsRateDecodingHPC_POP_MaxCorr.results.condition== cond);
                            rows = intersect(rows,cols);
                            mse_rate = lsRateDecodingHPC_POP_MaxCorr.results.mse_rate(rows);
                            mse_chance_rate = lsRateDecodingHPC_POP_MaxCorr.results.mse_chance_rate(rows);
                            fits_rate = lsRateDecodingHPC_POP_MaxCorr.results.fits(rows);
                            for trial = 1:length(fits_rate)
                               se_r(trial,:) = makeLength((fits_rate(trial).rate-fits_rate(trial).response').^2,1000); 
            %                    se_p(trial,:) = makeLength((fits_rate(trial).phase_all-fits_rate(trial).response).^2,1000);  
                               se_c(trial,:) = makeLength((fits_rate(trial).response(randperm(length(fits_rate(trial).response)))-fits_rate(trial).response).^2,1000);  
                            end
                            error_rate(count,:) = mean(se_r);
            %                 error_phase(count,:) = mean(se_p);
                            error_chance(count,:) = mean(se_c);
                            animal(count) = sum(double(sessionInfo.animal));
                            pvals(count,cond) = pval;
                            pvals_rate(count,cond) = pval_r;
                            
                            clear se_*
                            if isempty(sessionInfo.ca3)
                                region(count) = 1;
                            else
                                region(count) = 3;
                            end
                            figure(region(count))
                            subplot(2,2,1)
                            imagesc(error_rate(region==region(count),:))
                            subplot(2,2,2)
            %                 imagesc(error_phase(region==region(count),:))
                            subplot(2,2,3)
                            plot(mean(error_rate(region==region(count),:)),'r')
                            hold on
            %                 plot(mean(error_phase(region==region(count),:)),'g')
                            plot(mean(error_chance(region==region(count),:)),'k')
                            hold off
                            count=1+count;
%                         end
%                     end
                    end
            end
            end        
            end
            end
        end
    end
    
    pause(.01)
    cd D:\Dropbox\datasets\lsDataset
end