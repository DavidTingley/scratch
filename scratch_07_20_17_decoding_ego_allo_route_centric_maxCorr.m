% d  = dir('*201*');
% % ii=38;
% 
% for ii=1:68
%     cd(d(ii).name)
xml = LoadParameters;
load([xml.FileName '.behavior.mat'])
load([xml.FileName '.sessionInfo.mat'])
spikes = bz_GetSpikes;

%     if ~exist([spikes.sessionName '.referenceFramesMaxCorr.mat'])
if isfield(behavior.events.trials{1},'direction')
% lfp = bz_GetLFP(sessionInfo.thetaChans(2));
%     
conditions = unique(behavior.events.trialConditions);
nCells = length(spikes.times);
routeCentricSamplingRate = behavior.samplingRate;
smoothingRange = 1:50;

% % find a better way to get spike phase relationship...
% [firingMaps] = bz_firingMap1D(spikes,behavior,lfp,4);
load([spikes.sessionName '.firingMaps.cellinfo.mat'])
load([spikes.sessionName '.phaseMaps.cellinfo.mat'])
% % iterate through conditions and compile spike trains and spike-phase
% % trainse
routeCentric = [];
conditionCentric = [];
goalCentric = [];
alloCentric = [];
egoCentric = [];
for cond = conditions
    trials = find(behavior.events.trialConditions==cond);
    intervals = behavior.events.trialIntervals(trials,:);
    for t = 1:length(trials)
        trial = trials(t);
        if isfield(behavior.events.trials{trial},'direction')
        spk_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000)); % assumes intervals are in seconds, rounds to nearest millisecond
        phase_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000));
        for cell = 1:nCells
            if ~isempty(phaseMaps.phaseMaps{cond}{cell})
                f = find(phaseMaps.phaseMaps{cond}{cell}(:,2)==t);
                if ~isempty(f)
                for s=1:length(f)
                    phase_trains{cond}{t}(cell,ceil(phaseMaps.phaseMaps{cond}{cell}(f(s),5)*1000)) = ...
                        phaseMaps.phaseMaps{cond}{cell}(f(s),end);
                end
                end
            end
            sp = find(InIntervals(spikes.times{cell},intervals(t,:)));
            if ~isempty(sp)
                spks = Restrict(ceil((spikes.times{cell}(sp)-intervals(t,1))*1000+.0000001),[1 size(spk_trains{cond}{t},2)]);
                spk_trains{cond}{t}(cell,spks)=1;
            end
        end 
        
        nBins = size(spk_trains{cond}{t},2);
        nPos = length(behavior.events.trials{trial}.x);
        
        % goal centric is the distance to the reward for all trials.
        goalCentric = [goalCentric;makeLength(behavior.events.trials{trial}.mapping,201)'];
                % the -1 gaurantees the length to be longer than the above spk/phase trains

        % route centric, given rotations...
        if strcmp(behavior.events.trials{trial}.direction,'counter-clockwise')
            if strcmp(behavior.events.trials{trial}.type,'central alternation') 
                routeLoc = 1;
            elseif strcmp(behavior.events.trials{trial}.type,'wheel alternation')
                routeLoc = 2;
            elseif strcmp(behavior.events.trials{trial}.type,'linear')
                 routeLoc = 3;   
            end
        else
            if strcmp(behavior.events.trials{trial}.type,'central alternation') 
                routeLoc = 4;
            elseif strcmp(behavior.events.trials{trial}.type,'wheel alternation')
                routeLoc = 5;
            elseif strcmp(behavior.events.trials{trial}.type,'linear')
                routeLoc = 6;
            end
        end
        route = zeros(201,6);
        route(:,routeLoc) = makeLength(behavior.events.trials{trial}.mapping,201);
        routeCentric = [routeCentric; route];
        
        %
        c = zeros(201,length(conditions));
        c(:,cond) = makeLength(behavior.events.trials{trial}.mapping,201);
        conditionCentric = [conditionCentric; c];
        % and allocentric locations in space
        try
        angles = quat2eul([behavior.events.trials{trial}.orientation.x,...
            behavior.events.trials{trial}.orientation.y,...
            behavior.events.trials{trial}.orientation.z,...
            behavior.events.trials{trial}.orientation.w]);
        catch
        angles = [zeros(length(behavior.events.trials{trial}.orientation.yaw),1),...
            behavior.events.trials{trial}.orientation.yaw,...
            zeros(length(behavior.events.trials{trial}.orientation.yaw),1),...
            zeros(length(behavior.events.trials{trial}.orientation.yaw),1)];
        end
        try
        allo = [makeLength(behavior.events.trials{trial}.x,201)...
            ;makeLength(behavior.events.trials{trial}.y,201)...
            ;makeLength(behavior.events.trials{trial}.z,201)...
            ;acos(makeLength(cos(angles(:,1)),201))...
            ;acos(makeLength(cos(angles(:,2)),201))...
            ;acos(makeLength(cos(angles(:,3)),201))]';
        catch
        allo = [makeLength(behavior.events.trials{trial}.x,201)...
            ;makeLength(behavior.events.trials{trial}.y,201)...
            ;makeLength(zeros(length(behavior.events.trials{trial}.y),1),201)...
            ;acos(makeLength(cos(angles(:,1)),201))...
            ;acos(makeLength(cos(angles(:,2)),201))...
            ;acos(makeLength(cos(angles(:,3)),201))]';    
        end
        vels = [zeros(1,6); diff(allo)];
        acc =  [zeros(1,6); diff(diff(allo)); zeros(1,6)];
        egoCentric = [egoCentric; acc , vels];
        alloCentric = [alloCentric;allo]; clear a
        end
    end
end

% now set up neural data
warning off


for cell = 1:nCells
    c=1;
%     figure
for wind = smoothingRange
    maps = bz_firingMap1D(spikes,behavior,wind);
    [binnedPhaseMap] = bz_phaseMap2Bins(phaseMaps.phaseMaps,firingMaps.rateMaps,behavior);
%     for cell = 1:nCells 
    %smoothing
%         phase_trains_smooth=[];
%         rates_trains_smooth = [];
%         for cond = conditions
%             for t = 1:length(phase_trains{cond})
%             phase_trains_smooth=[phase_trains_smooth;...
%             circ_smoothTS(phase_trains{cond}{t}(cell,:),wind,'method','mean','exclude',0)];
%             rates_trains_smooth = [rates_trains_smooth; ...
%             smoothts(spk_trains{cond}{t}(cell,:),'b',wind)'];
%             end
%         end
            phase_trains_smooth=[];
            rates_trains_smooth = [];
            for cond = conditions
            for t = 1:size(binnedPhaseMap{cond},2)
                binnedPhaseMap{cond}(cell,t,isnan(binnedPhaseMap{cond}(cell,t,:)))=0;
                phase_trains_smooth=[phase_trains_smooth;...
                    circ_smoothTS(squeeze(binnedPhaseMap{cond}(cell,t,:)),wind,'method','mean','exclude',0)];
                rates_trains_smooth = [rates_trains_smooth; ...
                                       squeeze(maps.rateMaps_box{cond}(cell,t,:))];                   
            end 
            end
        
        % split data into train/test
        predictors = [alloCentric, routeCentric, conditionCentric, goalCentric, egoCentric];
                     % 6 6 2 1 12               
        for iter = 1:5
            r = randperm(length(rates_trains_smooth));
            pct = round(prctile(1:length(r),60));
            rates_train = rates_trains_smooth(r(1:pct));
            rates_test = rates_trains_smooth(r(pct+1:end));
            
            phase_train = phase_trains_smooth(r(1:pct));
            phase_test = phase_trains_smooth(r(pct+1:end));
            predictors_train = predictors(r(1:pct),:);
            predictors_test = predictors(r(pct+1:end),:);
            for p=1:size(predictors,2)
                if sum(rates_train)>100
            %% begin modelling rate
%             mse_chance_rate(c,iter) = mean((rates_test(randperm(length(rates_test)))-rates_test).^2);
            
%                 [b dev_p stats] = glmfit(predictors_train(:,p),rates_train,'normal');
%                 yfit = glmval(b,predictors_test(:,p),'identity');
%                 mse_rate(c,p,iter) = mean((yfit-rates_test).^2);     
                
            cl = max_correlation_coefficient_CL;
            cl = train(cl,predictors_train(:,p)',rates_train');
            yfit = test(cl,[predictors_test(:,p)']);
            mse_rate(c,p,iter) = mean((yfit-rates_test').^2);
     
            %% begin modelling phase
%             mse_chance_phase(c,iter) = mean((phase_test(randperm(length(phase_test)))-phase_test).^2);
%             mse_chance_phase_cos(c,iter) = mean((cos(phase_test(randperm(length(phase_test))))-cos(phase_test)).^2);
%             mse_chance_phase_sin(c,iter) = mean((sin(phase_test(randperm(length(phase_test))))-sin(phase_test)).^2);
%                 [b dev_p stats] = glmfit(predictors_train(:,p),phase_train,'normal');
%                 yfit = glmval(b,predictors_test(:,p),'identity');
%                 mse_phase(c,p,iter) = mean((yfit-phase_test).^2);   
%                 [b dev_p stats] = glmfit(predictors_train(:,p),cos(phase_train),'normal');
%                 yfit = glmval(b,predictors_test(:,p),'identity');
%                 mse_phase_cos(c,p,iter) = mean((yfit-cos(phase_test)).^2);   
%                 [b dev_p stats] = glmfit(predictors_train(:,p),sin(phase_train),'normal');
%                 yfit = glmval(b,predictors_test(:,p),'identity');
%                 mse_phase_sin(c,p,iter) = mean((yfit-sin(phase_test)).^2);   
                
                cl = max_correlation_coefficient_CL;
                cl = train(cl,predictors_train(:,p)',[(phase_train)]');
                yfit = test(cl,[predictors_test(:,p)']);
                mse_phase(c,p,iter) = mean((yfit-phase_test').^2);
                
                cl = max_correlation_coefficient_CL;
                cl = train(cl,predictors_train(:,p)',[cos(phase_train)]');
                yfit = test(cl,[predictors_test(:,p)']);
                mse_phase_cos(c,p,iter) = mean((yfit-cos(phase_test)').^2);
                
                cl = max_correlation_coefficient_CL;
                cl = train(cl,predictors_train(:,p)',[sin(phase_train)]');
                yfit = test(cl,[predictors_test(:,p)']);
                mse_phase_sin(c,p,iter) = mean((yfit-sin(phase_test)').^2);
%                  [b dev_p stats] = glmfit([phase_train cos(phase_train) sin(phase_train)],predictors_train(:,p),'normal');
%                 yfit = glmval(b,[phase_test cos(phase_test) sin(phase_test)],'identity');
%                 mse_phase(c,p,iter) = mean((yfit-rates_test).^2); 
                    
            else
                mse_rate(c,p,iter) = nan;
                mse_phase(c,p,iter) = nan;
                mse_phase_cos(c,p,iter) = nan;
                mse_phase_sin(c,p,iter) = nan;
                end
            end
        end
        
%         subplot(2,2,2)
%         for j=1:size(mse_rate,1)
%                 mse_norm(j,:,:) = zscore(mse_rate(j,:,:));
%         end
%         imagesc(1:p,smoothingRange(1:c),(squeeze(nanmean(mse_norm,3))))
%         title('rate')
%         clear mse_norm
%         
%         subplot(2,2,4)
%         for j=1:size(mse_phase_cos,1)
%                 mse_norm(j,:,:) = zscore(mse_phase_cos(j,:,:));
%         end
%         imagesc(1:p,smoothingRange(1:c),(squeeze(nanmean(mse_norm,3))))
%         clear mse_norm
%         title('phase')
%         hold off
%         
%         subplot(2,2,1)
%         % allo
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_rate(:,1:6,:),2),3)))),'.k')
%         hold on
%         %route
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_rate(:,7:10,:),2),3)))),'.r')
%         %cond
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_rate(:,11:20,:),2),3)))),'.b')
%         %goal
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_rate(:,21,:),2),3)))),'.g')
%         %ego
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_rate(:,22:end,:),2),3)))),'.m')
%         hold off
%         
%         subplot(2,2,3)
%         % allo
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,1:6,:),2),3))))./mean(mse_chance_phase,2),'.k')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_cos(:,1:6,:),2),3)))),'.k')
%         hold on
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_sin(:,1:6,:),2),3)))),'.k')
%         %route
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,7:10,:),2),3))))./mean(mse_chance_phase,2),'.r')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_cos(:,7:10,:),2),3)))),'.r')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_sin(:,7:10,:),2),3)))),'.r')
%         %cond
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,11:20,:),2),3))))./mean(mse_chance_phase,2),'.b')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_cos(:,11:20,:),2),3)))),'.b')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_sin(:,11:20,:),2),3)))),'.b')
%         %goal
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,21,:),2),3))))./mean(mse_chance_phase,2),'.g')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_cos(:,21,:),2),3)))),'.g')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_sin(:,21,:),2),3)))),'.g')
%         %ego
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,22:end,:),2),3))))./mean(mse_chance_phase,2),'.m')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_cos(:,22:end,:),2),3)))),'.m')
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_phase_sin(:,22:end,:),2),3)))),'.m')
%         hold off
%         
% %         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse,3))))
%         title(cell)
%         
%         pause(.05)
        c=c+1;
end
%     mse_all_phase{cell} = mse_phase;
    mse_all_phase_cos{cell} = mse_phase_cos;
    mse_all_phase_sin{cell} = mse_phase_sin;
    mse_all_phase{cell} = mse_phase;
    mse_all_rate{cell} = mse_rate;
%     mse_all_chance_phase_cos{cell} = mse_chance_phase_cos;
%     mse_all_chance_phase_sin{cell} = mse_chance_phase_sin;
%     mse_all_chance_rate{cell} = mse_chance_rate;
    clear mse mse_allo mse_cond mse_route mse_goal mse_chance_rate mse_chance_phase_cos mse_chance_phase_sin mse_rate mse_phase*
    save([spikes.sessionName '.referenceFramesMaxCorr.mat'],'mse*')
end
end
%     end
% cd /home/david/datasets/lsDataset
% end




