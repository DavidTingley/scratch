function [] = scratch_09_26_17_GLM_ls_hpc_populationDecoding_RATE_ONLY_from_ratemaps(RECORDING,septalCell)

cd([RECORDING])
plotting = 0;
saveMat = 0;
smoothingRange = [10 25 35 60 250 1000];

ls = bz_GetSpikes('region','ls');
spikes = bz_GetSpikes;

if septalCell <= length(ls.times)
   
    sessionInfo = bz_getSessionInfo;
load([sessionInfo.FileName '.behavior.mat'])
load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
conditions = unique(behavior.events.trialConditions);


for i = conditions
    for k=1:size(firingMaps.rateMaps{i},2)
%         septalResponse{i}{k} = phase_trains{i}{k}(septalCell,:);
        septalResponse{i}{k} = squeeze(firingMaps.rateMaps{i}(septalCell,k,:));
    end
end

for i= conditions
    for k=1:length(septalResponse{i})
        spk_trains{i}{k} = squeeze(firingMaps.rateMaps{i}(length(ls.times)+1:end,k,:));
        position{i}{k} = 1:201;
    end
end
nCells = length(spikes.times) - (length(ls.times));
%% set up MaxCorr decoder
lsRateDecodingHPC_POP_MaxCorr.results = table;

for cond = conditions
r = randperm(length(spk_trains{cond}));
for iter = 1:length(r)
r = [r(end) r(1:end-1)];r(end)
for wind = smoothingRange
        for cell = 1:nCells
            rates_trains_smooth = [];
            septalResponse_train = [];
            position_train = [];
            theta_train = [];
            for t = 1:length(spk_trains{cond})-1
                rates_trains_smooth = [rates_trains_smooth; smooth(spk_trains{cond}{r(t)}(cell,:),wind)*wind];
                septalResponse_train = [septalResponse_train; septalResponse{cond}{r(t)}];
                position_train = [position_train; position{cond}{r(t)}'];
%                 theta_train = [theta_train, theta{cond}{r(t)}];
            end
          
            rates_trains_smooth_train(cell,:) = rates_trains_smooth;
          
            rates_trains_smooth_test(cell,:) = [...
                smooth(spk_trains{cond}{r(end)}(cell,:),wind)*wind];
            
        end
        septalResponse_test = septalResponse{cond}{r(end)};
%         theta_test = theta{cond}{r(end)};
        
        %% rate coding model
        
        [res dev] = glmfit([rates_trains_smooth_train]',round(septalResponse_train),'normal');
        yfit_rate = glmval(res,[rates_trains_smooth_test;]','identity');
        
%         cl = max_correlation_coefficient_CL;
%         cl = train(cl,[rates_trains_smooth_train; theta_train],round(septalResponse_train));
        
%         for ts = 1:length(septalResponse_test)
%            yfit_rate(ts) = test(cl,[(rates_trains_smooth_test(:,ts));theta_test(ts)]); 
%         end
        struct.mse_rate = mean((yfit_rate'-septalResponse_test).^2);
        % chance rate
        rr = randperm(length(septalResponse_train));
        rrr = randperm(length(septalResponse_test));
        
        [res dev] = glmfit([rates_trains_smooth_train(:,rr);]',round(septalResponse_train),'normal');
        yfit_chance_rate = glmval(res,[rates_trains_smooth_test(:,rrr);]','identity');
        
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
save([RECORDING '/' sessionInfo.FileName '_rateMap_decoding_GLM_cell_' num2str(septalCell) '.mat'],'lsRateDecodingHPC_POP_MaxCorr')

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