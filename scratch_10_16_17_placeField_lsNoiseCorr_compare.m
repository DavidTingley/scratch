cd D:\Dropbox\datasets\lsDataset
clear all
count = 1;
d  = dir('*201*');
placeFieldPctThresh = 1; % 1 - percent diff b/w peak and trough to be defined as PF

 
%% compile data

for ii=1:length(d)
    cd(d(ii).name)
    if ~isempty(dir('*firingMaps*'))  & exist([d(ii).name '.noiseCorrs.mat']) & exist([d(ii).name '.placeFields.' num2str(placeFieldPctThresh,'%0.2i') '_pctThresh.mat']) & exist([d(ii).name '.positionDecodingGLM_binnedspace_box.cellinfo.mat'])
        disp(['working on ' d(ii).name ''])
        sessionInfo = bz_getSessionInfo;
        load([d(ii).name '.noiseCorrs.mat'],'noiseCorr')
        load([d(ii).name '.positionDecodingGLM_binnedspace_box.cellinfo.mat'])
        load([d(ii).name '.behavior.mat'],'behavior')
        load([d(ii).name '.placeFields.' num2str(placeFieldPctThresh,'%0.2i') '_pctThresh.mat'],'fields') % use PF defs with diff thresholds
        load([d(ii).name '.firingMaps.cellinfo.mat'],'firingMaps')
        spikes = bz_GetSpikes;
        for cell=1:length(spikes.times)
            if strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca3') | strcmp(spikes.region{cell},'ca1')
            if size(noiseCorr,2) >= cell
                for ls = 1:size(noiseCorr,1)
                    for cond = 1:size(noiseCorr,3)
                        % check that LS cell phase codes...
%                         rows = find(positionDecodingGLM_binnedspace_box.results{ls}.condition == cond);
%                         cols = find(positionDecodingGLM_binnedspace_box.results{ls}.tau == 50);
%                         rows = intersect(rows,cols);
%                         if ~isempty(rows)
%                         [a b] = kstest2(positionDecodingGLM_binnedspace_box.results{ls}.mse_phase_all(rows),positionDecodingGLM_binnedspace_box.results{ls}.mse_chance);
                        % check here for trial type (linear, central,
                        % wheel)
                        if strcmp(behavior.events.conditionType{cond},'central')% && b < .01
%                         fields{cond} = bz_getPlaceFields1D(firingMaps.rateMaps{cond},'minPeakRate',2,'percentThresh',.1);
        
                        if sum(sum(firingMaps.countMaps{cond}(ls,:,:))) > 10  % ls neuron has to fire this many spikes to bother 
                        if ~isempty(fields{cond}{cell}) & size(noiseCorr,2) >= cell % didnt fill out matrix with 0's
                           nc(count,:) = squeeze(noiseCorr(ls,cell,cond,:));  
                           meanRate(count,:) = squeeze(mean(firingMaps.rateMaps{cond}(cell,:,:)));
                           com(count) = fields{cond}{cell}{1}.COM;
                           animal(count) = sum(double(sessionInfo.animal));
                           if ~isempty(sessionInfo.ca3)
                              region(count) = 3;
                           else
                              region(count) = 1;
                           end
                           subplot(2,2,1);
                           plot(squeeze(noiseCorr(ls,cell,cond,:)));
                           subplot(2,2,2);
                           plot(meanRate(count,:))
                           subplot(2,2,3);
                           plot(squeeze(mean(firingMaps.rateMaps{cond}(ls,:,:))))
                           pause
                           count=1+count;
                        end
%                         end
                        end
                        end
                    end
                end
            end
            end
        end
    end
    cd D:\Dropbox\datasets\lsDataset
   
end


%% plotting
% remove nan/0/1 first?
nc(nc==1)=nan;
nc(nc==0)=nan;

for i=1:200
    if i < 21 
        ii = 21;
    elseif i > 180
        ii = 180;
    else
        ii=i;
    end
[a b] = find(abs(com-i)<20);
c1 = intersect(b,find(region==1));
c3 = intersect(b,find(region==3));
subplot(2,2,1)
errorbar(i,nanmean(nanmean(nc(c1,ii-20:ii+20))'),nanstd(nanmean(nc(c1,ii-20:ii+20))'),'.k')
hold on
errorbar(i,nanmean(nanmean(nc(c3,ii-20:ii+20))'),nanstd(nanmean(nc(c3,ii-20:ii+20))'),'.r')
ca1(i) = nanmean(nanmean(nc(c1,ii-20:ii+20))');
ca3(i) = nanmean(nanmean(nc(c3,ii-20:ii+20))');

subplot(2,2,2)
plot(i,nanmean(nanmean(nc(c1,ii-20:ii+20))')-nanmean(nanmean(nc(c3,ii-20:ii+20))'),'.g')
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



