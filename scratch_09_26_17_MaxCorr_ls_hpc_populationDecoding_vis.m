function [] = scratch_09_26_17_MaxCorr_ls_hpc_populationDecoding_vis(RECORDING,septalCell,COND)

cd([RECORDING])
plotting = 1;
saveMat = 0;
% smoothingRange = [2 4 6 10 20 25 40 60 120 250 500 1000 3000];
smoothingRange = [5 40 300];
ls = bz_GetSpikes('region','ls');
if septalCell <= length(ls.times)
    
spikes = bz_GetSpikes;
sessionInfo = bz_getSessionInfo;
load([sessionInfo.FileName '.behavior.mat'])
lfp = bz_GetLFP(sessionInfo.thetaChans(end));
% septalCell=2;

lsDecodingHPC_POP_MaxCorr.UID = spikes.UID;
lsDecodingHPC_POP_MaxCorr.region = spikes.region;
lsDecodingHPC_POP_MaxCorr.sessionName = spikes.sessionName; % session name

conditions = unique(behavior.events.trialConditions);
nCells = length(spikes.times);
positionSamplingRate = behavior.samplingRate;

% find a better way to get spike phase relationship...
[firingMaps] = bz_firingMap1D(spikes,behavior,lfp,5);
[b a] = butter(3,[4/626 12/625],'bandpass');
theta_phases = angle(hilbert(FiltFiltM(b,a,double(lfp.data))));
% iterate through conditions and compile spike trains and spike-phase
% trains
for cond = conditions
    trials = find(behavior.events.trialConditions==cond);
    intervals = behavior.events.trialIntervals(trials,:);
    for t = 1:length(trials)
        trial = trials(t);
        spk_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000)); % assumes intervals are in seconds, rounds to nearest millisecond
        phase_trains{cond}{t} = zeros(nCells,ceil((intervals(t,2)-intervals(t,1))*1000));
        
        for cell = 1:nCells
            if ~isempty(firingMaps.phaseMaps{cond}{cell})
                f = find(firingMaps.phaseMaps{cond}{cell}(:,2)==t);
                if ~isempty(f)
                for s=1:length(f)
                    phase_trains{cond}{t}(cell,ceil(firingMaps.phaseMaps{cond}{cell}(f(s),5)*1000)) = ...
                        firingMaps.phaseMaps{cond}{cell}(f(s),end);
                end
                end
            end
            sp = find(InIntervals(spikes.times{cell},intervals(t,:)));
            spk_trains{cond}{t}(cell,ceil((spikes.times{cell}(sp)-intervals(t,1))*1000+.00001))=1;
        end 
        
%         position{cond}{t} = interp1(1:length(behavior.events.trials{trial}.x)...
%             ,behavior.events.trials{trial}.mapping,1:positionSamplingRate/1000:length(...
%             behavior.events.trials{trial}.x));
        nBins = size(spk_trains{cond}{t},2);
        nPos = length(behavior.events.trials{trial}.x);
        position{cond}{t} = interp1(1:length(behavior.events.trials{trial}.x)...
            ,behavior.events.trials{trial}.mapping,1:(nPos-1)/nBins:length(...
            behavior.events.trials{trial}.x)); % the -1 gaurantees the length to be longer than the above spk/phase trains
        position{cond}{t} = position{cond}{t}(1:length(spk_trains{cond}{t}));
        take = FindInInterval(lfp.timestamps,intervals(t,:));
        theta{cond}{t} = makeLength(theta_phases(take(1):take(2)),length(position{cond}{t}));
    end
end

for i = conditions
    for k=1:length(position{i})
        septalResponse{i}{k} = spk_trains{i}{k}(septalCell,:);
        septalResponse_phase{i}{k} = phase_trains{i}{k}(septalCell,:);
    end
end

for i= conditions
    for k=1:length(septalResponse{i})
        spk_trains{i}{k} = spk_trains{i}{k}(length(ls.times)+1:end,:);
        phase_trains{i}{k} = phase_trains{i}{k}(length(ls.times)+1:end,:);
    end
end
nCells = nCells - (length(ls.times));
%% set up MaxCorr decoder
lsDecodingHPC_POP_MaxCorr.results = table;

for cond = COND
r = randperm(length(phase_trains{cond}));
for iter = 1:length(r)
r = circshift(r,1);r(end)
for wind = smoothingRange
        for cell = 1:nCells
            phase_trains_smooth=[];
            rates_trains_smooth = [];
            septalResponse_train = [];septalResponse_phase_train = [];
            position_train = [];
            theta_train = [];
            for t = 1:length(phase_trains{cond})-1
                phase_trains_smooth=[phase_trains_smooth; circ_smoothTS(phase_trains{cond}{r(t)}(cell,:),wind,'method','mean','exclude',0)];
                rates_trains_smooth = [rates_trains_smooth; smooth(spk_trains{cond}{r(t)}(cell,:),wind)*wind];
                septalResponse_train = [septalResponse_train; (septalResponse{cond}{r(t)}')];
                septalResponse_phase_train = [septalResponse_phase_train; (septalResponse_phase{cond}{r(t)}')];
%                 septalResponse_train = [septalResponse_train; smooth(septalResponse{cond}{r(t)}',wind)];
%                 septalResponse_train = [septalResponse_train; circ_smoothTS(septalResponse{cond}{r(t)}',wind,'method','mean','exclude',0)];
                position_train = [position_train; position{cond}{r(t)}'];
                theta_train = [theta_train, theta{cond}{r(t)}];
            end
            phase_trains_smooth_train(cell,:) = phase_trains_smooth;
            rates_trains_smooth_train(cell,:) = rates_trains_smooth;
            phase_trains_smooth_test(cell,:)=[...
                circ_smoothTS(phase_trains{cond}{r(end)}(cell,:),wind,'method','mean','exclude',0)];
            rates_trains_smooth_test(cell,:) = [...
                smooth(spk_trains{cond}{r(end)}(cell,:),wind)*wind];
            
        end
        
        septalResponse_test = (septalResponse{cond}{r(end)});
%         septalResponse_test = smooth(septalResponse{cond}{r(end)},wind)';
%         septalResponse_test_phase = circ_smoothTS(septalResponse_phase{cond}{r(end)},wind,'method','mean','exclude',0)';
        septalResponse_test_phase = septalResponse_phase{cond}{r(end)};
        
        theta_test = theta{cond}{r(end)};
        
        %% rate coding model
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[rates_trains_smooth_train; theta_train],(septalResponse_train));
        
        for ts = 1:length(septalResponse_test)
           yfit_rate(ts) = test(cl,[(rates_trains_smooth_test(:,ts));theta_test(ts)]); 
        end
        struct.mse_rate = circ_mean(abs(yfit_rate-septalResponse_test)');
        % chance rate
        rr = randperm(length(theta_train));
        rrr = randperm(length(theta_test));
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[rates_trains_smooth_train(:,rr); theta_train],(septalResponse_train));
        
        for ts = 1:length(septalResponse_test)
           yfit_chance_rate(ts) = test(cl,[(rates_trains_smooth_test(:,rrr(ts)));theta_test(ts)]); 
        end
        struct.mse_chance_rate = circ_mean(abs(yfit_chance_rate-septalResponse_test)');
        
        %% phase coding model
        
        % discretize phase_trains here...
        phase_trains_smooth_train_cos = cos(phase_trains_smooth_train);
        phase_trains_smooth_train_sin = sin(phase_trains_smooth_train);
        phase_trains_smooth_test_cos = cos(phase_trains_smooth_test);
        phase_trains_smooth_test_sin = sin(phase_trains_smooth_test);
        
%         phase_trains_smooth_train(phase_trains_smooth_train==0)=nan;
%         phase_trains_smooth_test(phase_trains_smooth_test==0)=nan;
%         phase_trains_smooth_train = discretize(phase_trains_smooth_train,-pi:.1:pi);
%         phase_trains_smooth_test = discretize(phase_trains_smooth_test,-pi:.1:pi);
%         phase_trains_smooth_train(isnan(phase_trains_smooth_train))=0;
%         phase_trains_smooth_test(isnan(phase_trains_smooth_test))=0;
%         
%         phase_trains_smooth_train_cos(phase_trains_smooth_train_cos==0)=nan;
%         phase_trains_smooth_test_cos(phase_trains_smooth_test_cos==0)=nan;
%         phase_trains_smooth_train_cos = discretize(phase_trains_smooth_train_cos,-1:.1:1);
%         phase_trains_smooth_test_cos = discretize(phase_trains_smooth_test_cos,-1:.1:1);
%         phase_trains_smooth_train_cos(isnan(phase_trains_smooth_train_cos))=0;
%         phase_trains_smooth_test_cos(isnan(phase_trains_smooth_test_cos))=0;
%         
%         phase_trains_smooth_train_sin(phase_trains_smooth_train_sin==0)=nan;
%         phase_trains_smooth_test_sin(phase_trains_smooth_test_sin==0)=nan;
%         phase_trains_smooth_train_sin = discretize(phase_trains_smooth_train_sin,-1:.1:1);
%         phase_trains_smooth_test_sin = discretize(phase_trains_smooth_test_sin,-1:.1:1);
%         phase_trains_smooth_train_sin(isnan(phase_trains_smooth_train_sin))=0;
%         phase_trains_smooth_test_sin(isnan(phase_trains_smooth_test_sin))=0;
        
        % non-transformed circular decoding
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[phase_trains_smooth_train; theta_train],(septalResponse_train));
        
        for ts = 1:length(septalResponse_test)
           yfit_circ(ts) = test(cl,[phase_trains_smooth_test(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase = circ_mean(abs(yfit_circ-septalResponse_test)');
        
        % cos and sin transformed circular decoding
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[phase_trains_smooth_train_cos; theta_train],(septalResponse_train));
        
        for ts = 1:length(septalResponse_test)
           yfit_cos(ts) = test(cl,[phase_trains_smooth_test_cos(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase_cos = circ_mean(abs(yfit_cos-septalResponse_test)');
        
        cl = max_correlation_coefficient_CL;
        cl = train(cl,[phase_trains_smooth_train_sin; theta_train],(septalResponse_train));
        
        for ts = 1:length(septalResponse_test)
           yfit_sin(ts) = test(cl,[phase_trains_smooth_test_sin(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase_sin = circ_mean(abs(yfit_sin-septalResponse_test)');
        
        % all phase models in one...
        cl = max_correlation_coefficient_CL;
        all_phase_train = [phase_trains_smooth_train_cos;phase_trains_smooth_train_sin];
        all_phase_test =  [phase_trains_smooth_test_cos;phase_trains_smooth_test_sin];
        cl = train(cl,[all_phase_train; theta_train],(septalResponse_train));
        
        for ts = 1:length(septalResponse_test)
           yfit_circ_all(ts) = test(cl,[all_phase_test(:,ts);theta_test(ts)]); 
        end
        struct.mse_phase_all = circ_mean(abs(yfit_circ_all-septalResponse_test)');
        
        % chance phase
        cl = max_correlation_coefficient_CL;
        all_phase_train = [phase_trains_smooth_train_cos;phase_trains_smooth_train_sin];
        all_phase_test =  [phase_trains_smooth_test_cos;phase_trains_smooth_test_sin];
        cl = train(cl,[all_phase_train(:,rr); theta_train],(septalResponse_train));
        
        for ts = 1:length(septalResponse_test)
           yfit_chance(ts) = test(cl,[all_phase_test(:,rrr(ts));theta_test(ts)]); 
        end
        struct.mse_chance = circ_mean(abs(yfit_chance-septalResponse_test)');
        fits = cell2struct({yfit_chance,yfit_rate,yfit_circ_all},{'chance','rate','phase_all'},2);
        %% put data into struct/table
        struct.tau = wind;
        struct.condition = cond;
        struct.iter = iter;
        struct.trialOrder = r;
        struct.fits = fits;
        struct.response = septalResponse_test;
        lsDecodingHPC_POP_MaxCorr.results = [lsDecodingHPC_POP_MaxCorr.results;struct2table(struct)];
        if plotting
            clf
            subplot(3,2,1)
             t_rate = varfun(@mean,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_rate',...
            'GroupingVariables',{'tau','condition'});
            t_phase = varfun(@mean,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_phase_all',...
            'GroupingVariables',{'tau','condition'});
            t_chance = varfun(@mean,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_chance',...
            'GroupingVariables',{'tau','condition'});
            t_chance_s = varfun(@std,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_chance',...
            'GroupingVariables',{'tau','condition'});
            t_chance_rate = varfun(@mean,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_chance_rate',...
            'GroupingVariables',{'tau','condition'});
            t_chance_s_rate = varfun(@std,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_chance_rate',...
            'GroupingVariables',{'tau','condition'});
            t_rate_s = varfun(@std,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_rate',...
            'GroupingVariables',{'tau','condition'});
            t_phase_s = varfun(@std,lsDecodingHPC_POP_MaxCorr.results,'InputVariables','mse_phase_all',...
            'GroupingVariables',{'tau','condition'});
            tab = join(join(t_rate,t_phase),t_chance);
            tab_s = join(join(t_rate_s,t_phase_s),t_chance_s);
            rows = find(tab.condition==cond);

            title('MaxCorr decoding of pos, r-rate, g-phase')
%             iterRows = lsDecodingHPC_POP_MaxCorr.iter == iter;
%             plot(lsDecodingHPC_POP_MaxCorr.results.tau,...
%                 lsDecodingHPC_POP_MaxCorr.results.mse_rate,'r')
%             hold on
%             plot(lsDecodingHPC_POP_MaxCorr.results.tau,...
%                 lsDecodingHPC_POP_MaxCorr.results.mse_phase_all,'g')
%             hold off
            boundedline(tab.tau,t_chance_rate.mean_mse_chance_rate(rows),t_chance_s_rate.std_mse_chance_rate(rows),'k')
            boundedline(tab.tau,tab.mean_mse_chance(rows),tab_s.std_mse_chance(rows),'k')
            boundedline(tab.tau,tab.mean_mse_phase_all(rows),tab_s.std_mse_phase_all(rows),'g')
            boundedline(tab.tau,tab.mean_mse_rate(rows),tab_s.std_mse_rate(rows),'r')
            set(gca,'xscale','log')
            subplot(3,2,2)
            plot(lsDecodingHPC_POP_MaxCorr.results.tau,...
                lsDecodingHPC_POP_MaxCorr.results.mse_phase_cos,'g')
            subplot(3,2,3)
            plot(lsDecodingHPC_POP_MaxCorr.results.tau,...
                lsDecodingHPC_POP_MaxCorr.results.mse_phase_sin,'g')
            subplot(3,2,4)
            plot(lsDecodingHPC_POP_MaxCorr.results.tau,...
                lsDecodingHPC_POP_MaxCorr.results.mse_phase,'g')
            subplot(3,2,5)
            
            for i=1:201
                f = find((abs(position{cond}{r(end)}-i)<10));
                ff = find((abs(position{cond}{r(end)}-i)<10) & septalResponse_test_phase~=0);
                if ~isempty(f)
%                     d = discretize(yfit_rate(f),2).*theta_test(f);
%                     mean_rate(iter,i) = circ_mean([d(d~=0)]');
%                     d = discretize(yfit_circ_all(f),2).*theta_test(f);
%                     mean_phase(iter,i) = circ_mean([d(d~=0)]');
                    d = (yfit_circ_all(f));
                    mean_phase(iter,i) = circ_mean([d(d~=0)]');
                    d = (yfit_rate(f));
                    mean_rate(iter,i) = circ_mean([d(d~=0)]');
                    
                    mean_actual(iter,i) = circ_mean([septalResponse_test_phase(ff)]');
                else
                    mean_rate(iter,i)=nan;
                    mean_phase(iter,i)  = nan;
                    mean_actual(iter,i) = nan;
                end
            end
            plot(unwrap(circ_mean(mean_actual)),'k');
            hold on
            plot(unwrap(circ_mean(mean_phase)),'g');
            plot(unwrap(circ_mean(mean_rate)),'r');
            plot(unwrap(circ_mean(mean_actual))+2*pi,'k');
            plot(unwrap(circ_mean(mean_phase))+2*pi,'g');
            plot(unwrap(circ_mean(mean_rate))+2*pi,'r');
            axis([0 201 -pi pi*3])
            hold off
            
            subplot(3,2,6)
            
            errorbar(2,mean(lsDecodingHPC_POP_MaxCorr.results.mse_rate),std(lsDecodingHPC_POP_MaxCorr.results.mse_rate),'r')
            hold on
            errorbar(3,mean(lsDecodingHPC_POP_MaxCorr.results.mse_phase_all),std(lsDecodingHPC_POP_MaxCorr.results.mse_phase_all),'g')
            errorbar(1,mean(lsDecodingHPC_POP_MaxCorr.results.mse_chance),std(lsDecodingHPC_POP_MaxCorr.results.mse_chance),'k')
            [a pc] = kstest2(lsDecodingHPC_POP_MaxCorr.results.mse_phase_all,lsDecodingHPC_POP_MaxCorr.results.mse_chance);
            [a rc] = kstest2(lsDecodingHPC_POP_MaxCorr.results.mse_rate,lsDecodingHPC_POP_MaxCorr.results.mse_chance);
            [a rp] = kstest2(lsDecodingHPC_POP_MaxCorr.results.mse_rate,lsDecodingHPC_POP_MaxCorr.results.mse_phase);
            errorbar(4,mean(lsDecodingHPC_POP_MaxCorr.results.mse_chance_rate),std(lsDecodingHPC_POP_MaxCorr.results.mse_chance_rate),'k')
            axis([0 5 0 9])
            title([ num2str(pc) '   '  num2str(rc) '    ' num2str(rp)])
%             sigstar({[1,2], [1,3], [2,3]},[rc pc rp])
            hold off
            pause(.001)
        end
        clear *train *test yfit*
        disp(['finished with window: ' num2str(wind) ' out of ' num2str(smoothingRange(end)) ' total'])
end

end
save([RECORDING '/' sessionInfo.FileName '_cell_' num2str(septalCell) '.mat'],'lsDecodingHPC_POP_MaxCorr')
% 
% errorbar(2,mean(lsDecodingHPC_POP_MaxCorr.results.mse_rate),std(lsDecodingHPC_POP_MaxCorr.results.mse_rate),'r')
% hold on
% errorbar(3,mean(lsDecodingHPC_POP_MaxCorr.results.mse_phase_all),std(lsDecodingHPC_POP_MaxCorr.results.mse_phase_all),'g')
% errorbar(1,mean(lsDecodingHPC_POP_MaxCorr.results.mse_chance),std(lsDecodingHPC_POP_MaxCorr.results.mse_chance),'k')
% [a pc] = kstest2(lsDecodingHPC_POP_MaxCorr.results.mse_phase_all,lsDecodingHPC_POP_MaxCorr.results.mse_chance)
% [a rc] = kstest2(lsDecodingHPC_POP_MaxCorr.results.mse_rate,lsDecodingHPC_POP_MaxCorr.results.mse_chance)
% [a rp] = kstest2(lsDecodingHPC_POP_MaxCorr.results.mse_rate,lsDecodingHPC_POP_MaxCorr.results.mse_phase)
% errorbar(4,mean(lsDecodingHPC_POP_MaxCorr.results.mse_chance_rate),std(lsDecodingHPC_POP_MaxCorr.results.mse_chance_rate),'k')
% axis([0 5 0 9])
% title([ num2str(pc) '   '  num2str(rc) '    ' num2str(rp)])
% sigstar({[1,2], [1,3], [2,3]},[rc pc rp])
end
end