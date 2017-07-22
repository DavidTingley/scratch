% d  = dir('*201*');
% % ii=38;
% 
% for ii=1:length(d)
xml = LoadParameters;
load([xml.FileName '.behavior.mat'])
load([xml.FileName '.sessionInfo.mat'])
spikes = bz_GetSpikes;
lfp = bz_GetLFP(sessionInfo.thetaChans(2));
%     
conditions = unique(behavior.events.trialConditions);
nCells = length(spikes.times);
routeCentricSamplingRate = behavior.samplingRate;
smoothingRange = 1:300:4000;
% 
% % find a better way to get spike phase relationship...
[rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,4);
% 
% % iterate through conditions and compile spike trains and spike-phase
% % trains
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
        spk_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000)); % assumes intervals are in seconds, rounds to nearest millisecond
        phase_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000));
        for cell = 1:nCells
            if ~isempty(phaseMap{cond}{cell})
                f = find(phaseMap{cond}{cell}(:,2)==t);
                if ~isempty(f)
                for s=1:length(f)
                    phase_trains{cond}{t}(cell,ceil(phaseMap{cond}{cell}(f(s),5)*1000)) = ...
                        phaseMap{cond}{cell}(f(s),end);
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
        goalCentric = [goalCentric;makeLength(behavior.events.trials{trial}.mapping,length(spk_trains{cond}{t}))'];
                % the -1 gaurantees the length to be longer than the above spk/phase trains

        % route centric, given rotations...
        if strcmp(behavior.events.trials{trial}.direction,'counter-clockwise')
            if strcmp(behavior.events.trials{trial}.type,'central alternation')
                routeLoc = 1;
            else
                routeLoc = 2;
            end
        else
            if strcmp(behavior.events.trials{trial}.type,'central alternation')
                routeLoc = 3;
            else
                routeLoc = 4;
            end
        end
        route = zeros(length(spk_trains{cond}{t}),4);
        route(:,routeLoc) = makeLength(behavior.events.trials{trial}.mapping,length(spk_trains{cond}{t}));
        routeCentric = [routeCentric; route];
        
        %
        c = zeros(length(spk_trains{cond}{t}),length(conditions));
        c(:,cond) = makeLength(behavior.events.trials{trial}.mapping,length(spk_trains{cond}{t}));
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
        allo = [makeLength(behavior.events.trials{trial}.x,length(spk_trains{cond}{t}))...
            ;makeLength(behavior.events.trials{trial}.y,length(spk_trains{cond}{t}))...
            ;makeLength(behavior.events.trials{trial}.z,length(spk_trains{cond}{t}))...
            ;acos(makeLength(cos(angles(:,1)),length(spk_trains{cond}{t})))...
            ;acos(makeLength(cos(angles(:,2)),length(spk_trains{cond}{t})))...
            ;acos(makeLength(cos(angles(:,3)),length(spk_trains{cond}{t})))]';
        catch
        allo = [makeLength(behavior.events.trials{trial}.x,length(spk_trains{cond}{t}))...
            ;makeLength(behavior.events.trials{trial}.y,length(spk_trains{cond}{t}))...
            ;makeLength(zeros(length(behavior.events.trials{trial}.y),1),length(spk_trains{cond}{t}))...
            ;acos(makeLength(cos(angles(:,1)),length(spk_trains{cond}{t})))...
            ;acos(makeLength(cos(angles(:,2)),length(spk_trains{cond}{t})))...
            ;acos(makeLength(cos(angles(:,3)),length(spk_trains{cond}{t})))]';    
        end
        vels = [zeros(1,6); diff(allo)];
        acc =  [zeros(1,6); diff(diff(allo)); zeros(1,6)];
        egoCentric = [egoCentric; acc , vels];
        alloCentric = [alloCentric;allo]; clear a
    end
end

% now set up neural data
warning off


for cell =1:nCells
    c=1;
%     figure
for wind = smoothingRange
%     for cell = 1:nCells 
    %smoothing
        phase_trains_smooth=[];
        rates_trains_smooth = [];
        for cond = conditions
            for t = 1:length(phase_trains{cond})
            phase_trains_smooth=[phase_trains_smooth;...
            circ_smoothTS(phase_trains{cond}{t}(cell,:),wind,'method','mean','exclude',0)];
            rates_trains_smooth = [rates_trains_smooth; ...
            smoothts(spk_trains{cond}{t}(cell,:),'b',wind)'];
            end
        end
        
        % split data into train/test
        predictors = [alloCentric, routeCentric, conditionCentric, goalCentric, egoCentric];
        
        for iter = 1:5
            r = randperm(length(rates_trains_smooth));
            pct = round(prctile(1:length(r),60));
            rates_train = rates_trains_smooth(r(1:pct));
            rates_test = rates_trains_smooth(r(pct+1:end));
            
            phase_train = phase_trains_smooth(r(1:pct));
            phase_test = phase_trains_smooth(r(pct+1:end));
            predictors_train = predictors(r(1:pct),:);
            predictors_test = predictors(r(pct+1:end),:);

            %% begin modelling rate
            % allo
%             [b dev_allo stats] = glmfit([predictors_train(:,1:6)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,1:6),'identity');
%             mse_allo(c,iter) = mean((yfit-rates_test).^2); 
%             % route
%             [b dev_route stats] = glmfit([predictors_train(:,7:10)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,7:10),'identity');
%             mse_route(c,iter) = mean((yfit-rates_test).^2); 
%             % condition
%             [b dev_cond stats] = glmfit([predictors_train(:,11:20)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,11:20),'identity');
%             mse_cond(c,iter) = mean((yfit-rates_test).^2); 
%             % goal
%             [b dev_goal stats] = glmfit([predictors_train(:,21)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,21),'identity');
%             mse_goal(c,iter) = mean((yfit-rates_test).^2); 
%             %ego
%             [b dev_goal stats] = glmfit([predictors_train(:,22:end)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,22:end),'identity');
%             mse_ego(c,iter) = mean((yfit-rates_test).^2); 
            
            mse_chance_rate(c,iter) = mean((rates_test(randperm(length(rates_test)))-rates_test).^2);
            for p=1:size(predictors,2)
                [b dev_p stats] = glmfit(predictors_train(:,p),rates_train,'normal');
                yfit = glmval(b,predictors_test(:,p),'identity');
                mse_rate(c,p,iter) = mean((yfit-rates_test).^2);     
            end
            %% begin modelling phase
            % allo
%             [b dev_allo stats] = glmfit([predictors_train(:,1:6)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,1:6),'identity');
%             mse_allo(c,iter) = mean((yfit-rates_test).^2); 
%             % route
%             [b dev_route stats] = glmfit([predictors_train(:,7:10)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,7:10),'identity');
%             mse_route(c,iter) = mean((yfit-rates_test).^2); 
%             % condition
%             [b dev_cond stats] = glmfit([predictors_train(:,11:20)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,11:20),'identity');
%             mse_cond(c,iter) = mean((yfit-rates_test).^2); 
%             % goal
%             [b dev_goal stats] = glmfit([predictors_train(:,21)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,21),'identity');
%             mse_goal(c,iter) = mean((yfit-rates_test).^2); 
%             %ego
%             [b dev_goal stats] = glmfit([predictors_train(:,22:end)],rates_train,'normal');
%             yfit = glmval(b,predictors_test(:,22:end),'identity');
%             mse_ego(c,iter) = mean((yfit-rates_test).^2); 
            
            mse_chance_phase(c,iter) = mean((phase_test(randperm(length(phase_test)))-phase_test).^2);
            for p=1:size(predictors,2)
                [b dev_p stats] = glmfit(predictors_train(:,p),phase_train,'normal');
                yfit = glmval(b,predictors_test(:,p),'identity');
                mse_phase(c,p,iter) = mean((yfit-phase_test).^2);     
            end            
            
        end
        
%         subplot(2,2,1)
%         plot(smoothingRange(1:c),mean(mse_allo,2),'.k')
%         hold on
%         plot(smoothingRange(1:c),mean(mse_route,2),'.r')
%         plot(smoothingRange(1:c),mean(mse_cond,2),'.b')
%         plot(smoothingRange(1:c),mean(mse_goal,2),'.g')
%         set(gca,'yscale','log')
%         hold off
%         subplot(2,2,2)
%         for j=1:size(mse_rate,1)
%                 mse_norm(j,:,:) = zscore(mse_rate(j,:,:));
%         end
%         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse_norm,3))))
%         title('rate')
%         clear mse_norm
%         subplot(2,2,4)
%         for j=1:size(mse_phase,1)
%                 mse_norm(j,:,:) = zscore(mse_phase(j,:,:));
%         end
%         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse_norm,3))))
%         clear mse_norm
%         title('phase')
%         hold off
%         subplot(2,2,1)
%         % allo
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,1:6,:),2),3))))./mean(mse_chance_rate,2),'.k')
%         hold on
%         %route
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,7:10,:),2),3))))./mean(mse_chance_rate,2),'.r')
%         %cond
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,11:20,:),2),3))))./mean(mse_chance_rate,2),'.b')
%         %goal
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,21,:),2),3))))./mean(mse_chance_rate,2),'.g')
%         %ego
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,22:end,:),2),3))))./mean(mse_chance_rate,2),'.m')
%         hold off
%         
%         subplot(2,2,3)
%         % allo
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,1:6,:),2),3))))./mean(mse_chance_phase,2),'.k')
%         hold on
%         %route
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,7:10,:),2),3))))./mean(mse_chance_phase,2),'.r')
%         %cond
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,11:20,:),2),3))))./mean(mse_chance_phase,2),'.b')
%         %goal
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,21,:),2),3))))./mean(mse_chance_phase,2),'.g')
%         %ego
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,22:end,:),2),3))))./mean(mse_chance_phase,2),'.m')
%         hold off
%         
% %         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse,3))))
%         title(cell)
% %         subplot(2,2,2)
% %         [b dev stats] = glmfit([alloCentric, conditionCentric, routeCentric, goalCentric],sin(phase_trains_smooth),'normal');
% %         plot(wind,dev,'.k')
% %         hold on
% %         [b dev stats] = glmfit([alloCentric, goalCentric],phase_trains_smooth,'normal');
% %         plot(wind,dev,'.r')
% %         [b dev stats] = glmfit([routeCentric, goalCentric],phase_trains_smooth,'normal');
% %         plot(wind,dev,'.g')
% %         [b dev stats] = glmfit([alloCentric, routeCentric],phase_trains_smooth,'normal');
% %         plot(wind,dev,'.m')
% % %         hold off
% %         set(gca,'yscale','log')
%         pause(.05)
        c=c+1;
end
    mse_all_phase{cell} = mse_phase;
    mse_all_rate{cell} = mse_rate;
    mse_chance_phase{cell} = mse_chance_phase;
    mse_chance_rate{cell} = mse_chance_rate;
    clear mse mse_allo mse_cond mse_route mse_goal mse_chance
end

save([spikes.sessionName '.referenceFrames.mat'],'mse*')
% end




