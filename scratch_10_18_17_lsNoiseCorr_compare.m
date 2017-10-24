cd /home/david/Dropbox/datasets/lsDataset
clear all
count = 1;
cellCount = 1;
d  = dir('*201*');
placeFieldPctThresh = 20; % 1 - percent diff b/w peak and trough to be defined as PF
celltotal = 1;
 
%% compile data

for ii=1:length(d)
    cd(d(ii).name)
    if ~isempty(dir('*firingMaps*'))  & exist([d(ii).name '.noiseCorrs.mat']) & exist([d(ii).name '.placeFields.' num2str(placeFieldPctThresh,'%0.2i') '_pctThresh.mat'])
        disp(['working on ' d(ii).name ''])
        sessionInfo = bz_getSessionInfo;
        load([d(ii).name '.noiseCorrs.mat'],'noiseCorr')
        load([d(ii).name '.behavior.mat'],'behavior')
        load([d(ii).name '.positionDecodingGLM_binnedspace_box.cellinfo.mat'])
%         load([d(ii).name '.placeFields.' num2str(placeFieldPctThresh,'%0.2i') '_pctThresh.mat'],'fields') % use PF defs with diff thresholds
        load([d(ii).name '.firingMaps.cellinfo.mat'],'firingMaps')
        spikes = bz_GetSpikes;
        for cell=1:length(spikes.times)
            if strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca3') | strcmp(spikes.region{cell},'ca1')
            if size(noiseCorr,2) >= cell
                for ls = 1:size(noiseCorr,1)
                    for cond = 1:size(noiseCorr,3)
                        % check if LS cell phase codes..
                        rows = find(positionDecodingGLM_binnedspace_box.results{ls}.condition == cond);
                        cols = find(positionDecodingGLM_binnedspace_box.results{ls}.tau == 50);
                        rows = intersect(rows,cols);
                        if ~isempty(rows)
                        [a b] = kstest2(positionDecodingGLM_binnedspace_box.results{ls}.mse_phase_all(rows),positionDecodingGLM_binnedspace_box.results{ls}.mse_chance);
                        % check here for trial type (linear, central,
                        % wheel)
                        pvals(cellCount,cond) = b;
                        if strcmp(behavior.events.conditionType{cond},'wheel') && b < .01
%                         fields{cond} = bz_getPlaceFields1D(firingMaps.rateMaps{cond},'minPeakRate',2,'percentThresh',.1);
        
                        if sum(sum(firingMaps.countMaps{cond}(ls,:,:))) > 20  % ls neuron has to fire this many spikes to bother 
%                         if ~isempty(fields{cond}{cell}) & size(noiseCorr,2) >= cell % didnt fill out matrix with 0's
                           nc(count,:) = squeeze(noiseCorr(ls,cell,cond,:)); 
                           meanRate(count,:) = squeeze(mean(firingMaps.rateMaps{cond}(cell,:,:)));
%                            com(count) = fields{cond}{cell}{1}.COM;
                           signalCorr(count) = corr(meanRate(count,:)',squeeze(mean(firingMaps.rateMaps{cond}(ls,:,:))));
                           animal(count) = sum(double(sessionInfo.animal));
                           precess(count) = 1;
                           
                           if ~isempty(sessionInfo.ca3)
                              region(count) = 3;
                           else
                              region(count) = 1;
                           end
                           count=1+count;
%                         end
                        end
                        end
                        end
                    end
                end
            end
            end
        end
    end
    cd /home/david/Dropbox/datasets/lsDataset
   
end


%% plotting
% remove nan/0/1 first?
nc(nc>.999)=nan;
nc(nc==0)=nan;

for i=1:200
    if i < 41 
        ii = 41;
    elseif i > 160
        ii = 160;
    else
        ii=i;
    end
% [a b] = find(abs(com-i)<20);
c1 = find(region==1);
c3 = find(region==3);
subplot(2,2,1)
% errorbar(i,nanmedian(nanmedian(nc(c1,i))'),nanstd(nanmedian(nc(c1,i))'),'.k')
plot(i,nanmedian(nanmedian(nc(c1,i))'),'.k')
hold on
errorbar(i,nanmedian(nanmedian(nc(c3,i))'),nanstd(nanmedian(nc(c3,i))'),'.r')
plot(i,nanmedian(nanmedian(nc(c3,i))'),'.r')
ca1(i) = nanmedian(nanmedian(nc(c1,i))');
ca3(i) = nanmedian(nanmedian(nc(c3,i))');

subplot(2,2,2)
plot(i,nanmedian(nanmedian(nc(c1,i))')-nanmedian(nanmedian(nc(c3,i))'),'.g')
hold on
end
for i=1:size(meanRate,1)
    r(i,:) = minmax_norm(meanRate(i,:));
end
subplot(2,2,3)
[a b o] =sort_cells(r(region==1,:),r(region==1,:),1);
imagesc(a)
title(size(a,1))
subplot(2,2,4)
[a b o] =sort_cells(r(region==3,:),r(region==3,:),1);
imagesc(a)
title(size(a,1))



