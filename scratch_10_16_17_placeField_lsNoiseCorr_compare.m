cd /home/david/datasets/lsDataset
clear all
count = 1;
d  = dir('*201*');
placeFieldPctThresh = 10; % 1 - percent diff b/w peak and trough to be defined as PF

 
%% compile data

for ii=1:length(d)
    cd(d(ii).name)
    if ~isempty(dir('*firingMaps*')) & exist([d(ii).name '.placeFields.' num2str(placeFieldPctThresh,'%0.2i') '_pctThresh.mat']) & exist([d(ii).name '.noiseCorrs.mat'])
        disp(['working on ' d(ii).name ''])
        sessionInfo = bz_getSessionInfo;
        load([d(ii).name '.noiseCorrs.mat'],'noiseCorr')
        load([d(ii).name '.placeFields.' num2str(placeFieldPctThresh,'%0.2i') '_pctThresh.mat'],'fields') % use PF defs with diff thresholds
        load([d(ii).name '.firingMaps.cellinfo.mat'],'firingMaps')
        spikes = bz_GetSpikes;
        for cell=1:length(spikes.times)
            if strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca3') | strcmp(spikes.region{cell},'ca1')
                for ls = 1:size(noiseCorr,1)
                    for cond = 1:size(noiseCorr,3)
                        % check here for trial type (linear, central,
                        % wheel)

                        if sum(sum(firingMaps.countMaps{cond}(ls,:,:))) > 20  % ls neuron has to fire this many spikes to bother 
                        if ~isempty(fields{cond}{cell}) & size(noiseCorr,2) >= cell % didnt fill out matrix with 0's
                           nc(count,:) = squeeze(noiseCorr(ls,cell,cond,:)); 
                           meanRate(count,:) = squeeze(mean(firingMaps.rateMaps{cond}(cell,:,:)));
                           com(count) = fields{cond}{cell}{1}.COM;
                           if ~isempty(sessionInfo.ca3)
                              region(count) = 3;
                           else
                              region(count) = 1;
                           end
                           count=1+count;
                        end
                        end
                    end
                end
            end
        end
    end
    cd /home/david/datasets/lsDataset
   
end


%% plotting
for i=11:192
[a b] = find(abs(com-i)<20);
c1 = intersect(b,find(region==1));
c3 = intersect(b,find(region==3));
subplot(2,2,1)
plot(i,nanmean(nanmean(nc(c1,i-10:i+10))'),'.k')
hold on
plot(i,nanmean(nanmean(nc(c3,i-10:i+10))'),'.r')
subplot(2,2,2)
plot(i,nanmean(nanmean(nc(c1,i-10:i+10))')-nanmean(nanmean(nc(c3,i-10:i+10))'),'.g')
hold on
end
