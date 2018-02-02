d = dir('*201*');
count=1;
for rec = 1:length(d)
    cd(d(rec).name)
    if exist([d(rec).name '.firingMaps.cellinfo.mat'])
    load([d(rec).name '.phaseMaps.cellinfo.mat'])
    load([d(rec).name '.firingMaps.cellinfo.mat'])
    load([d(rec).name '.placeFields.20_pctThresh.mat'])
    load([d(rec).name '.behavior.mat'])
    
%% get velocity
for trial = 1:length(behavior.events.trials)
    vel{trial} = makeLength(smooth(abs(diff(behavior.events.trials{trial}.x))+abs(diff(behavior.events.trials{trial}.y)),50),200);
%     vel{trial} = mean(smooth(abs(diff(behavior.events.trials{trial}.x))+abs(diff(behavior.events.trials{trial}.y)),40),200);
end

%% iterate through conditions
    for cell=1:length(phaseMaps.UID)
    if ~strcmp(phaseMaps.region{cell},'ls')
    for cond = 1:length(phaseMaps.phaseMaps)
%         f = find(behavior.events.trialConditions==cond);
%         for t = 1:length(f)
%             trial = f(t);
        for field = 1:length(fields{cond}{cell})
        if ~isempty(phaseMaps.phaseMaps{cond}{cell}) & ~isempty(fields{cond}{cell}) & length(fields{cond}{cell})==1
        spkCount = 1;
        for spk = 1:size(phaseMaps.phaseMaps{cond}{cell},1)
            if phaseMaps.phaseMaps{cond}{cell}(spk,1) >= fields{cond}{cell}{field}.start && phaseMaps.phaseMaps{cond}{cell}(spk,1) <= fields{cond}{cell}{field}.stop
            trial = phaseMaps.phaseMaps{cond}{cell}(spk,2);
            bin = phaseMaps.phaseMaps{cond}{cell}(spk,1);
            dat{count}(spkCount,1) = vel{trial}(bin);
            dat{count}(spkCount,2) = (phaseMaps.phaseMaps{cond}{cell}(spk,end-1));%firingMaps.rateMaps{cond}(cell,trial,bin); %
            dat{count}(spkCount,3) = (phaseMaps.phaseMaps{cond}{cell}(spk,end));
            spkCount = 1 + spkCount;
            end
        end
        
        if size(dat{count},1) > 10
        phaseCorr(count) = circ_corrcl(dat{count}(:,3),dat{count}(:,1));
        rateCorr(count) = corr(dat{count}(:,2),dat{count}(:,1));
        for iter = 1:100
           r = randperm(size(dat{count},1));
           phaseCorr_shuffle(count,iter) = circ_corrcl(dat{count}(r,3),dat{count}(:,1));
           rateCorr_shuffle(count,iter) = corr(dat{count}(r,2),dat{count}(:,1));
        end
        else
            phaseCorr(count) = nan;
            rateCorr(count) = nan;
            phaseCorr_shuffle(count,:) = nan(100,1);
            rateCorr_shuffle(count,:) = nan(100,1);
        end
%         if t == 1
        count = 1+count;
%         end
        end
        end
    end
    end
    end
    else 
        warning(['not using: ' d(rec).name])
    end
    cd /home/david/datasets/lsDataset
end