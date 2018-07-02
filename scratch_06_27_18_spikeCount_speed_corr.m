
cd /home/david/datasets/lsDataset/
d = dir('*201*')
c=1;
v=[];
for f = 1:length(d)
    cd(d(f).name)
spikes = bz_GetSpikes('noprompt',true);
if ~isempty(spikes)
load([d(f).name '.firingMaps.cellinfo.mat'])
load([d(f).name '.behavior.mat'])
load([d(f).name '.placeFields.05_pctThresh.mat'])
    for i=1:length(spikes.times)
        if strcmp(spikes.region{i},'hpc') || strcmp(spikes.region{i},'ca1') || strcmp(spikes.region{i},'ca3')
        for j=1:length(firingMaps.rateMaps)
            if ~isempty(fields{j}{i})
            for t = 1:size(firingMaps.countMaps{j},2)
                if ~isempty(fields{j}{i})
                    start = fields{j}{i}{1}.start;
                    stop = fields{j}{i}{1}.stop;
                    nSpks{j}(i,t) = sum(firingMaps.countMaps{j}(i,t,start:stop));
                    rate{j}(i,t) = mean(firingMaps.rateMaps{j}(i,t,start:stop));
                    fieldIDX = find(ismember(behavior.events.trials{t}.mapping,start:stop));
                    if strcmp(behavior.units,'m')
                        velocity = (abs(diff(behavior.events.trials{t}.x))+abs(diff(behavior.events.trials{t}.y)))*100/(1/120);
                    else
                        velocity = (abs(diff(behavior.events.trials{t}.x))+abs(diff(behavior.events.trials{t}.y)))./10/(1/120);
                    end
                    vel{j}(i,t) = nanmean(velocity(fieldIDX)); %clear velocity
                else
                    vel{j}(i,t) = nan;
                    nSpks{j}(i,t) = nan;  
                    rate{j}(i,t) = nan;
                end
            end  
%             scatter(vel{j}(i,:),nSpks{j}(i,:),'.k')
            
            [corrs{j}(i) pval{j}(i)] = corr(vel{j}(i,:)',nSpks{j}(i,:)');
            [corrs_rate{j}(i) pval_rate{j}(i)] = corr(vel{j}(i,:)',rate{j}(i,:)');
            cc(c) = corrs{j}(i);
            pp(c) = pval{j}(i);
            rr(c) = corrs_rate{j}(i);
            pp_rate(c) = pval_rate{j}(i);
            v = [v,vel{j}(i,:)];
            c=1+c;    
            
            end
        end
        end
    end
    if exist('vel')
        velo{f} = vel;clear vel;
        numSpikes{f} = nSpks; clear nSpks;
        meanRates{f} = rate; clear rate;
        correlations{f} = corrs; clear corrs
        pvalues{f} = pval; clear pval
    end
end
    cd /home/david/datasets/lsDataset/
end