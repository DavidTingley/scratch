% try
    xml = LoadParameters;
     load([xml.FileName '.behavior.mat'])
    load([xml.FileName '.sessionInfo.mat'])
    lfp = bz_GetLFP(sessionInfo.thetaChans(end));%,'intervals',[behavior.timestamps(1) behavior.timestamps(end)]);
    [b a] = butter(4,[6/(lfp.samplingRate/2) 10/(lfp.samplingRate/2)],'bandpass');
    phases = angle(hilbert(FiltFiltM(b,a,double(lfp.data(:,1)))));
    
    spikes = bz_GetSpikes;
    
    pairs =[];
    for i=1:length((spikes.times))
       for j = i:length((spikes.times)) 
          if i ~= j & spikes.shankID(i) ~= spikes.shankID(j)
              if strcmp(spikes.region{i},'ls') && strcmp(spikes.region{j},'hpc')
                
                   pairs = [pairs; i j]; 
                   
              elseif strcmp(spikes.region{i},'hpc') && strcmp(spikes.region{j},'ls')
                   
                    pairs = [pairs; i j]; 
                    
              end
          end
       end
    end
  
    for c = 1:length(unique(behavior.events.trialConditions))
        spktrains{c} = [];
        spk_phase_trains{c} = [];
        phasetrains{c} =[];
        phasetrains_sin{c} = [];
        phasetrains_cos{c} =[];
        coords{c} = [];
        velocities{c} = [];
    end

    for c = 1:size(behavior.events.trialIntervals,1)
          trialLength = ceil(1000*(behavior.events.trialIntervals(c,2)-...
              behavior.events.trialIntervals(c,1)));

          train = single(zeros(length(spikes.times),trialLength));
          ptrain = single(zeros(length(spikes.times),trialLength));

          for s = 1:length(spikes.times)
                f = find(spikes.times{s}>behavior.events.trialIntervals(c,1));
                ff = find(spikes.times{s}<behavior.events.trialIntervals(c,2));
                fff = intersect(f,ff);
                times = ceil(1000*(spikes.times{s}(fff)-behavior.events.trialIntervals(c,1)));
                train(s,times) = 1;
                ptrain(s,times) = phases(round(spikes.times{s}(fff)*1250));
          end
          spktrains{behavior.events.trialConditions(c)} = ...
              [spktrains{behavior.events.trialConditions(c)},train];
          spk_phase_trains{behavior.events.trialConditions(c)} = ...
              [spk_phase_trains{behavior.events.trialConditions(c)},ptrain];
          [a start] =min(abs(lfp.timestamps-behavior.events.trialIntervals(c,1)));
          [a stop] =min(abs(lfp.timestamps-behavior.events.trialIntervals(c,2)));
          p = makeLength(phases(start:stop),length(train));
          ps = makeLength(sin(phases(start:stop)),length(train));
          pc = makeLength(cos(phases(start:stop)),length(train));

          phasetrains{behavior.events.trialConditions(c)} = ...
              [phasetrains{behavior.events.trialConditions(c)},p];
          phasetrains_sin{behavior.events.trialConditions(c)} = ...
              [phasetrains_sin{behavior.events.trialConditions(c)},ps];
          phasetrains_cos{behavior.events.trialConditions(c)} = ...
              [phasetrains_cos{behavior.events.trialConditions(c)},pc];

          cc = makeLength(behavior.events.trials{c}.mapping,length(train));
          coords{behavior.events.trialConditions(c)} = ...
              [coords{behavior.events.trialConditions(c)},cc];

          vel = abs(diff(behavior.events.trials{c}.x)) + abs(diff(behavior.events.trials{c}.y));
          v = makeLength(vel,length(train));
          velocities{behavior.events.trialConditions(c)} = ...
              [velocities{behavior.events.trialConditions(c)},v];
        % call GLM here
    end
        
 %% get phase coding stuff
 clear mark
for i=1:length(spikes.times)
    if strcmp(spikes.region{i},'ls')
        mark(i) = 1;
    else
        mark(i) = 0;
    end
end
ls = find(mark==1);


load([sessionInfo.FileName '.positionDecodingGLM_binnedspace_box.cellinfo.mat'])
load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
positionDecodingGLM = positionDecodingGLM_binnedspace_box;
conditions = length(unique(behavior.events.trialConditions));
nBins = round(length(behavior.events.map{1}.x));
for cell =1:length(ls)
    t_rate = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_rate',...
    'GroupingVariables',{'tau','condition'});
    t_phase = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_phase_all',...
    'GroupingVariables',{'tau','condition'});
    t_chance = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_chance',...
    'GroupingVariables',{'tau','condition'});
    tab = join(join(t_rate,t_phase),t_chance);
    t_phase_pval = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_phase_all_pval',...
    'GroupingVariables',{'tau','condition'});
    t_rate_pval = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_rate_pval',...
    'GroupingVariables',{'tau','condition'});
    pvals = join(t_phase_pval,t_rate_pval);
    for cond = 1:conditions
    if sum(behavior.events.trialConditions==cond) >= 12 %%%%%%%%%%%%%%%%%%%%%%%%%%
    nTrials = size(firingMaps.rateMaps{cond},2);
    %% information theory stuff here
    %% carry on..
    rows = find(tab.condition==cond);
    %                if sqrt(tab.mean_mse_phase_all(rows(1)))./nBins < .3
    [a b] =min(tab.mean_mse_phase_all(rows));
    [aa bb] =min(tab.mean_mse_rate(rows));
    first500ms = find(ismember(tab.tau(rows),1:nBins));
    %                 min_mse_rate = (min(tab.mean_mse_rate(rows(first500ms)))./mean(tab.mean_mse_chance(rows(first500ms))));
    %                 min_mse_phase_all = (min(tab.mean_mse_phase_all(rows(first500ms)))./mean(tab.mean_mse_chance(rows(first500ms))));
    min_mse_rate = (min(tab.mean_mse_rate(rows(first500ms)))./mean(tab.mean_mse_chance(rows)));
    min_mse_phase_all = (min(tab.mean_mse_phase_all(rows(first500ms)))./mean(tab.mean_mse_chance(rows)));
    min_mse_chance = (min(tab.mean_mse_chance(rows)))./mean(tab.mean_mse_chance(rows));
    NMSE_rates{cell}(cond) = min_mse_rate;
    NMSE_phases{cell}(cond) = min_mse_phase_all;
    NMSE_chance{cell}(cond) = min_mse_chance;
    end
    end
end
        
        
for i=1:length(spikes.times)
    if strcmp(spikes.region{i},'hpc')
        mark(i) = 1;
    else
        mark(i) = 0;
    end
end
f = find(mark==1);

for cond = 1:length(unique(behavior.events.trialConditions))
peers{cond} = spktrains{cond}(f,:);
peers_phase{cond} = spk_phase_trains{cond}(f,:);
    for iter = 1:10
       for ts = 1:size(peers_phase{cond},2)
           r = randperm(length(f));
           peers_phase_shuffle{cond}{iter}(:,ts) = peers_phase{cond}(r,ts);
       end
    end
end
clear *smooth
i=35;
for cond = 1:length(unique(behavior.events.trialConditions))
    for k = 1:size(peers{cond},1)
    peers_smooth{cond}(k,:) = smooth(peers{cond}(k,:),i);
    peers_phase_smooth{cond}(k,:) = circ_smoothTS(peers_phase{cond}(k,:),i,'method','mean','exclude',0);
    for iter = 1:10
        peers_phase_smooth_shuffle{cond}{iter}(k,:) = circ_smoothTS(peers_phase_shuffle{cond}{iter}(k,:),i,'method','mean','exclude',0);
    end
    end
end

 
 %% start the tuning models

for cell=1:length(ls)
    for cond = 1:length(unique(behavior.events.trialConditions))
        if sum(behavior.events.trialConditions==cond) >= 12 %%%%%%%%%%%%%%%%%%%%%%%%%%
        response = spktrains{cond}(cell,:);
        clear  dev dev_phase
        [results dev] = glmfit([peers_smooth{cond}',zscore(phasetrains{cond}'),zscore(phasetrains_cos{cond}'),zscore(phasetrains_sin{cond}'),zscore(round(coords{cond}))'],response,'normal');
        [results_phase dev_phase] = glmfit([cos(peers_phase_smooth{cond}'),zscore(phasetrains{cond}'),zscore(phasetrains_cos{cond}'),zscore(phasetrains_sin{cond}'),zscore(round(coords{cond}))'],response,'normal');
        for iter = 1:10
            [results_phase_shuffle{iter} dev_phase_shuffle(iter)] = glmfit([cos(peers_phase_smooth_shuffle{cond}{iter}'),zscore(phasetrains{cond}'),zscore(phasetrains_cos{cond}'),zscore(phasetrains_sin{cond}'),zscore(round(coords{cond}))'],response,'normal');
        end
        subplot(2,2,1)
        hold on
        plot(NMSE_phases{cell}(cond),dev-dev_phase,'.k')
        set(gca,'yscale','log')
        axis([0.3 1 .0000002 1])
        ylabel('HPC rate-phase tuning')
        xlabel('phase coding strength')
        subplot(2,2,2)
        plot(NMSE_rates{cell}(cond)-NMSE_phases{cell}(cond),dev-dev_phase,'.k'); hold on
        set(gca,'yscale','log')
        axis([-.5 .5 .0000002 1])
        subplot(2,2,3)
        plot(NMSE_rates{cell}(cond),dev-dev_phase,'.k'); hold on
        set(gca,'yscale','log')
        axis([0.3 1 .0000002 1])
        subplot(2,2,4)
        plot(dev_phase,mean(dev_phase_shuffle)-dev_phase,'.k')
        set(gca,'xscale','log')
        hold on
        pause(.1)
        
        devs_shuffle{cell}(cond,:) = dev_phase_shuffle;
        devs_rate{cell}(cond) = dev;
        devs_phase{cell}(cond) = dev_phase;
        end
    end
end  

save([sessionInfo.FileName '.ls_phase_rate_tuning.mat'],'devs*','NMSE*')

%         
%         
%         
% catch
%     end
%         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

