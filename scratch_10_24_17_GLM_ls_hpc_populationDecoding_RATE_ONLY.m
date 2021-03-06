function [] = scratch_09_26_17_GLM_ls_hpc_populationDecoding_RATE_ONLY(RECORDING,septalCell)

cd([RECORDING])
plotting = 0;
saveMat = 0;
smoothingRange = [10 25 35 60 250 1000];

ls = bz_GetSpikes('region','ls');
if septalCell <= length(ls.times)
    
spikes = bz_GetSpikes;
sessionInfo = bz_getSessionInfo;
load([sessionInfo.FileName '.behavior.mat'])
if isempty(sessionInfo.ca3)
    lfp = bz_GetLFP(sessionInfo.ca1);
else
    lfp = bz_GetLFP(sessionInfo.ca3);
end
% septalCell=2;

lsRateDecodingHPC_POP_MaxCorr.UID = spikes.UID;
lsRateDecodingHPC_POP_MaxCorr.region = spikes.region;
lsRateDecodingHPC_POP_MaxCorr.sessionName = spikes.sessionName; % session name

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
%         septalResponse{i}{k} = phase_trains{i}{k}(septalCell,:);
        septalResponse{i}{k} = spk_trains{i}{k}(septalCell,:);
    end
end

for i= conditions
    for k=1:length(septalResponse{i})
        spk_trains{i}{k} = spk_trains{i}{k}(length(ls.times)+1:end,:);
        phase_trains{i}{k} = phase_trains{i}{k}(length(ls.times)+1:end,:);
    end
end
nCells = nCells - (length(ls.times)+1);
%% set up MaxCorr decoder
lsRateDecodingHPC_POP_MaxCorr.results = table;

for cond = conditions
r = randperm(length(phase_trains{cond}));
for iter = 1:length(r)
r = [r(end) r(1:end-1)];r(end)
for wind = smoothingRange
        for cell = 1:nCells
            rates_trains_smooth = [];
            septalResponse_train = [];
            position_train = [];
            theta_train = [];
            for t = 1:length(phase_trains{cond})-1
                rates_trains_smooth = [rates_trains_smooth; smooth(spk_trains{cond}{r(t)}(cell,:),wind)*wind];
                septalResponse_train = [septalResponse_train; septalResponse{cond}{r(t)}'];
                position_train = [position_train; position{cond}{r(t)}'];
                theta_train = [theta_train, theta{cond}{r(t)}];
            end
          
            rates_trains_smooth_train(cell,:) = rates_trains_smooth;
          
            rates_trains_smooth_test(cell,:) = [...
                smooth(spk_trains{cond}{r(end)}(cell,:),wind)*wind];
            
        end
        septalResponse_test = septalResponse{cond}{r(end)};
        theta_test = theta{cond}{r(end)};
        
        %% rate coding model
        
        [res dev] = glmfit([rates_trains_smooth_train; cos(theta_train); sin(theta_train)]',round(septalResponse_train),'normal');
        yfit_rate = glmval(res,[rates_trains_smooth_test; cos(theta_test); sin(theta_test)]','identity');
        
%         cl = max_correlation_coefficient_CL;
%         cl = train(cl,[rates_trains_smooth_train; theta_train],round(septalResponse_train));
        
%         for ts = 1:length(septalResponse_test)
%            yfit_rate(ts) = test(cl,[(rates_trains_smooth_test(:,ts));theta_test(ts)]); 
%         end
        struct.mse_rate = mean((yfit_rate'-septalResponse_test).^2);
        % chance rate
        rr = randperm(length(theta_train));
        rrr = randperm(length(theta_test));
        
        [res dev] = glmfit([rates_trains_smooth_train(:,rr); cos(theta_train); sin(theta_train)]',round(septalResponse_train),'normal');
        yfit_chance_rate = glmval(res,[rates_trains_smooth_test(:,rrr); cos(theta_test); sin(theta_test)]','identity');
        
%         cl = max_correlation_coefficient_CL;
%         cl = train(cl,[rates_trains_smooth_train(:,rr); theta_train],round(septalResponse_train));
%         
%         for ts = 1:length(septalResponse_test)
%            yfit_chance_rate(ts) = test(cl,[(rates_trains_smooth_test(:,rrr(ts)));theta_test(ts)]); 
%         end
        struct.mse_chance_rate = mean((yfit_chance_rate'-septalResponse_test).^2);
        
     
        
        %% put data into struct/table
        fits = cell2struct({yfit_rate,septalResponse_test,yfit_chance_rate},{'rate','response','chance'},2);
        struct.tau = wind;
        struct.condition = cond;
        struct.iter = iter;
        struct.testTrial = r(end);
        struct.fits = fits;
        struct.septalCell = septalCell;
        lsRateDecodingHPC_POP_MaxCorr.results = [lsRateDecodingHPC_POP_MaxCorr.results;struct2table(struct)];
       
        clear *train *test yfit*
        disp(['finished with window: ' num2str(wind) ' out of ' num2str(smoothingRange(end)) ' total'])
end

end
save([RECORDING '/' sessionInfo.FileName '_rate_decoding_GLM_cell_' num2str(septalCell) '.mat'],'lsRateDecodingHPC_POP_MaxCorr')

% errorbar(2,mean(lsRateDecodingHPC_POP_MaxCorr.results.mse_rate),std(lsRateDecodingHPC_POP_MaxCorr.results.mse_rate),'r')
% hold on
% errorbar(3,mean(lsRateDecodingHPC_POP_MaxCorr.results.mse_phase_all),std(lsRateDecodingHPC_POP_MaxCorr.results.mse_phase_all),'g')
% errorbar(1,mean(lsRateDecodingHPC_POP_MaxCorr.results.mse_chance),std(lsRateDecodingHPC_POP_MaxCorr.results.mse_chance),'k')
% [a pc] = kstest2(lsRateDecodingHPC_POP_MaxCorr.results.mse_phase_all,lsRateDecodingHPC_POP_MaxCorr.results.mse_chance)
% [a rc] = kstest2(lsRateDecodingHPC_POP_MaxCorr.results.mse_rate,lsRateDecodingHPC_POP_MaxCorr.results.mse_chance)
% [a rp] = kstest2(lsRateDecodingHPC_POP_MaxCorr.results.mse_rate,lsRateDecodingHPC_POP_MaxCorr.results.mse_phase)
% errorbar(4,mean(lsRateDecodingHPC_POP_MaxCorr.results.mse_chance_rate),std(lsRateDecodingHPC_POP_MaxCorr.results.mse_chance_rate),'k')
% axis([0 5 0 9])
% title([ num2str(pc) '   '  num2str(rc) '    ' num2str(rp)])
% sigstar({[1,2], [1,3], [2,3]},[pc rc rp])
end
end