function [] =scratch_07_11_17_spatialBinnedBayes_population_ls(ensembleSize)
% if ~exist([sessionInfo.FileName '.posDecodeBayes.cellinfo.mat'])
warning off
sessionInfo = bz_getSessionInfo;
load([sessionInfo.FileName '.behavior.mat'])

% lfp = bz_GetLFP(sessionInfo.thetaChans(2));cd
nBins = length(behavior.events.map{1}.x);
% 

% % %    % set up phase coding data
% [firingMaps] = bz_firingMap1D(spikes,behavior,lfp,lfp4);
% if ~exist([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
    spikes = bz_GetSpikes('region','ls');
    r = randperm(length(spikes.times));
    ensemble = r(1:ensembleSize); 
    spikes = bz_GetSpikes('UID',spikes.UID(ensemble));clear r
    if ~exist(['/ifs/data/buzsakilab/popDecodeBayes/' sessionInfo.FileName '.posDecodeBayes.' num2str(sort(spikes.UID),'%0.3d') '.popinfo.mat'])
    load([sessionInfo.FileName '.behavior.mat'])
    if ~isempty(sessionInfo.ca1)
    lfp = bz_GetLFP(sessionInfo.ca1);
    elseif ~isempty(sessionInfo.ca3)
    lfp = bz_GetLFP(sessionInfo.ca3);
    else
    lfp = bz_GetLFP(sessionInfo.thetaChans(end));
    end
    if ~isempty(spikes)
    [firingMaps] = bz_firingMap1D(spikes,behavior,5,'saveMat',false);
    [phaseMaps] = bz_phaseMap1D(spikes,behavior,lfp,5,'saveMat',false);
    end
% else
%     load([sessionInfo.FileName '.firingMaps.cellinfo.mat']);
%     load([sessionInfo.FileName '.phaseMaps.cellinfo.mat']);
% end
rateMap =firingMaps.rateMaps;
countMap =firingMaps.countMaps;
occuMap =firingMaps.occupancy;
phaseMap=phaseMaps.phaseMaps;
[binnedPhaseMap] = bz_phaseMap2Bins(phaseMap,rateMap,behavior);
    
posDecodeBayes.region = spikes.region;
posDecodeBayes.sessionName = spikes.sessionName;
posDecodeBayes.UID = spikes.UID;
for cell=1:length(spikes.times)
    posDecodeBayes.results{cell} = table;
end

for smoothing = 1:round(nBins./2)
    disp(['smoothing by: ' num2str(smoothing) ' bins']);
    for cond = 1:length(unique(behavior.events.trialConditions))
        if size(binnedPhaseMap{cond},2) >= 5
%         figure(cond)
        % smooth data..
        maps = bz_firingMap1D(spikes,behavior,smoothing);
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
%                 rates_trains_smooth = [rates_trains_smooth; ...
%                                        smoothts(squeeze(countMap{cond}(cell,t,:))','b',smoothing)'];
                rates_trains_smooth = [rates_trains_smooth; ...
                                       squeeze(maps.rateMaps_box{cond}(cell,t,:))];   
              
                pos = [pos,[1:nBins]];                   
            end 
            r(cell,:) = rates_trains_smooth;
            p(cell,:) = phase_trains_smooth;
            p_cos(cell,:) = cos(p(cell,:));
            p_sin(cell,:) = sin(p(cell,:));
        end
        
            r=r';
            p=p';
            p_cos=p_cos';
            p_sin=p_sin'; 
            
            count = 1;
            for iter = 1:10
            rr = randperm(length(r));
            pct = round(prctile(1:length(r),60));
            
            r_train = discretize(r(rr(1:pct),:),min(r(rr(1:pct),:)):(max(r(rr(1:pct),:))-min(r(rr(1:pct),:)))./20:max(r(rr(1:pct),:)));
            r_test = discretize(r(rr(pct+1:end),:),min(r(rr(pct+1:end),:)):(max(r(rr(pct+1:end),:))-min(r(rr(pct+1:end),:)))./20:max(r(rr(pct+1:end),:)));
            p_train = discretize(p(rr(1:pct),:),min(p(rr(1:pct),:)):(max(p(rr(1:pct),:))-min(p(rr(1:pct),:)))./20:max(p(rr(1:pct),:)));
            p_test = discretize(p(rr(pct+1:end),:),min(p(rr(pct+1:end),:)):(max(p(rr(pct+1:end),:))-min(p(rr(pct+1:end),:)))./20:max(p(rr(pct+1:end),:)));
            p_cos_train = discretize(p_cos(rr(1:pct),:),min(p_cos(rr(1:pct),:)):(max(p_cos(rr(1:pct),:))-min(p_cos(rr(1:pct),:)))./20:max(p_cos(rr(1:pct),:)));
            p_cos_test = discretize(p_cos(rr(pct+1:end),:),min(p_cos(rr(pct+1:end),:)):(max(p_cos(rr(pct+1:end),:))-min(p_cos(rr(pct+1:end),:)))./20:max(p_cos(rr(pct+1:end),:)));
            p_sin_train = discretize(p_sin(rr(1:pct),:),min(p_sin(rr(1:pct),:)):(max(p_sin(rr(1:pct),:))-min(p_sin(rr(1:pct),:)))./20:max(p_sin(rr(1:pct),:)));
            p_sin_test = discretize(p_sin(rr(pct+1:end),:),min(p_sin(rr(pct+1:end),:)):(max(p_sin(rr(pct+1:end),:))-min(p_sin(rr(pct+1:end),:)))./20:max(p_sin(rr(pct+1:end),:)));
            
            
            pos_train = pos(rr(1:pct));
            pos_test = pos(rr(pct+1:end));
            % rate 
            cl = poisson_naive_bayes_CL;
            cl = train(cl,[r_train'],pos_train');
            yfit = test(cl,[r_test']);
            struct.mse_rate = nanmean((pos_test'-yfit).^2);
%             struct.mse_rate_pval = stats.p';
            % phase all
            cl = poisson_naive_bayes_CL;
            cl = train(cl,[p_cos_train'; p_sin_train';p_train'],pos_train');
            yfit = test(cl,[p_cos_test'; p_sin_test';p_test']);
            struct.mse_phase_all = nanmean((pos_test'-yfit).^2);
%             struct.mse_phase_all_pval = stats.p';
            %phase
            cl = poisson_naive_bayes_CL;
            cl = train(cl,[p_train'],pos_train');
            yfit = test(cl,[p_test']);
            struct.mse_phase = nanmean((pos_test'-yfit).^2);
%             struct.mse_phase_pval = stats.p';
            %phase cos
            cl = poisson_naive_bayes_CL;
            cl = train(cl,[p_cos_train'],pos_train');
            yfit = test(cl,[p_cos_test']);
            struct.mse_phase_cos = nanmean((pos_test'-yfit).^2);
%             struct.mse_phase_cos_pval = stats.p';
            %phase sin
%             cl = poisson_naive_bayes_CL;
%             cl = train(cl,[p_sin_train'],pos_train');
%             yfit = test(cl,[p_sin_test']);
%             struct.mse_phase_sin = nanmedian((pos_test'-yfit).^2);
%             struct.mse_phase_sin_pval = stats.p';
            % all predictors
            cl = poisson_naive_bayes_CL;
            cl = train(cl,[p_cos_train'; p_sin_train';p_train';r_train'],pos_train');
            yfit = test(cl,[p_cos_test'; p_sin_test';p_test';r_test']);
            struct.mse_both  = nanmean((pos_test'-yfit).^2);
%             struct.mse_both_pval = stats.p';
            
            cl = poisson_naive_bayes_CL;
            cl = train(cl,[r_train'],pos_train');
            yfit = test(cl,[r_test(randperm(length(r_test)),:)']);
            struct.mse_chance_rate  = nanmean((pos_test'-yfit).^2);
            
            cl = poisson_naive_bayes_CL;
            cl = train(cl,[p_cos_train'; p_sin_train';p_train'],pos_train');
            yfit = test(cl,[p_cos_test(randperm(length(p_test)),:)'; p_sin_test(randperm(length(p_test)),:)';p_test(randperm(length(p_test)),:)']);
            struct.mse_chance_phase  = nanmean((pos_test'-yfit).^2);
            
            % store peripherals
            struct.iter = count;
%             struct.dfe = stats.dfe;
            struct.tau = smoothing;
            struct.condition = cond;
            posDecodeBayes.results{cell} = [posDecodeBayes.results{cell}; struct2table(struct)];
            
            count = 1+count;
            end
            clear r p p_cos p_sin
           
% %                 
% % %                 figure(cond)
%                 t_rate = varfun(@mean,posDecodeBayes.results{cell},'InputVariables','mse_rate',...
%                 'GroupingVariables',{'tau','condition'});
%                 t_phase_all = varfun(@mean,posDecodeBayes.results{cell},'InputVariables','mse_phase_all',...
%                 'GroupingVariables',{'tau','condition'});
%                 t_phase = varfun(@mean,posDecodeBayes.results{cell},'InputVariables','mse_phase',...
%                 'GroupingVariables',{'tau','condition'});
%                 t_phase_cos = varfun(@mean,posDecodeBayes.results{cell},'InputVariables','mse_phase_cos',...
%                 'GroupingVariables',{'tau','condition'});
%                 t_both = varfun(@mean,posDecodeBayes.results{cell},'InputVariables','mse_both',...
%                 'GroupingVariables',{'tau','condition'});
%                 t_chance_rate = varfun(@mean,posDecodeBayes.results{cell},'InputVariables','mse_chance_rate',...
%                 'GroupingVariables',{'tau','condition'});   
%                 t_chance_phase = varfun(@mean,posDecodeBayes.results{cell},'InputVariables','mse_chance_phase',...
%                 'GroupingVariables',{'tau','condition'});  
%                 t_std_rate = varfun(@std,posDecodeBayes.results{cell},'InputVariables','mse_chance_rate',...
%                 'GroupingVariables',{'tau','condition'});  
%                 t_std_phase = varfun(@std,posDecodeBayes.results{cell},'InputVariables','mse_chance_phase',...
%                 'GroupingVariables',{'tau','condition'});  
%                 tab = join(join(join(join(join(join(t_rate,t_phase_all),t_both),t_chance_rate),t_phase),t_phase_cos),t_chance_phase);
%                 rows = find(tab.condition==cond);
%                 subplot(2,2,1);
%                 plot(tab.tau(rows),tab.mean_mse_phase(rows),'.g')
%                 hold on
%                 plot(tab.tau(rows),tab.mean_mse_phase_cos(rows),'g')
%                 hold off
% %                 imagesc(phase_trains_smooth);
% %                 caxis([-pi pi])
%                 subplot(2,2,2);
%                 boundedline(tab.tau(rows),tab.mean_mse_chance_rate(rows),3*t_std_rate.std_mse_chance_rate(rows),'k')
%                 hold on
% %                 nans(smoothing) = sum(isnan(yfit));
%                 plot(tab.tau(rows),tab.mean_mse_rate(rows),'r')
%                 plot(tab.tau(rows),tab.mean_mse_both(rows),'k')
%                 plot(tab.tau(rows),tab.mean_mse_phase_all(rows),'g')
%                 subplot(2,2,3);
%                 boundedline(tab.tau(rows),tab.mean_mse_chance_phase(rows),3*t_std_phase.std_mse_chance_phase(rows),'k')
%                 hold on
% 
% %                 plot(posDecodeBayes.results{cell}.tau(rows),...
% %                     mean(posDecodeBayes.results{cell}.mse_rate(rows,:)'),'r')
% %                 hold on
% %                 plot(posDecodeBayes.results{cell}.tau(rows),...
% %                     mean(posDecodeBayes.results{cell}.mse_phase_all(rows,:)'),'g')
%                 plot(posDecodeBayes.results{cell}.tau(rows),...
%                     mean(posDecodeBayes.results{cell}.mse_both(rows,:)'),'k')
%                 hold off
% %                 subplot(2,2,3)
% %                 imagesc(rates_trains_smooth);
% %                 subplot(2,2,4)
%                 
%                 title([cell cond smoothing]')
% 
%                 pause(.1)
            end         
        disp(['done with condition: ' num2str(cond) ' of ' num2str(length(unique(behavior.events.trialConditions)))]);
    end
    posDecodeBayes.dateRun = date;  % this can take a very long time so lets save each loop...
    posDecodeBayes.ensemble = ensemble;
    cd('/ifs/data/buzsakilab/popDecodeBayes')
    save(['sessionInfo.FileName '.' num2str(sort(spikes.UID),'%0.3d') '.popinfo.mat'],'posDecodeBayes')
end
end
