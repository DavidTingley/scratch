idx = find(pval(:,1)<.05 & corrs(:,1) > .1);

for i=1:length(idx)
rec = recording(idx(i));
idx_hpc = find(strcmp(content(rec).region{1},'hpc'));
if isempty(idx_hpc)
idx_hpc = find(strcmp(content(rec).region{1},'ca3'));
end


predictors = [nansum(content(rec).nSpikes{1}(idx_hpc,:))' ... %               [(nanmean(content(rec).PF{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0)))./nanmean(content(rec).nSpikes{1}(idx_hpc,:)~=0)]'... % % of spikes from PF%               [[nanmean(content(rec).rewardContent{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)))] ./ nanmean(content(rec).nSpikes{1}(idx_hpc,:))]'...
              content(rec).interRipInterval{1}',...
              nanmean(content(rec).cellLoc_wav{1}(idx_hpc,:).*content(rec).nSpikes{1}(idx_hpc,:))'...
              content(rec).condition{1}...
              content(rec).location{1}...
              content(rec).SleepState{1}...
              content(rec).duration{1}'...
              max(bay{rec})'...
              rankOrder_bayes{rec}'...
              rankOrder{rec}'...
              max(seqs{rec}')'];%...


actual = content(rec).nSpikes{1}(cellID(idx(i)),:);

keep = find(~isnan((predictors(:,9)))); 
keep2 = keep;
% keep = 1:size(predictors,1);

if ~isempty(unique(actual(keep)))
clf
subplot(2,2,1)
% scatter(actual(keep),predictors(keep,1),'.k')
keep = find(~isnan((predictors(:,2)))); 
boxplot(predictors(keep,2),actual(keep),'color','k')
hold on
plot(unique(actual(keep))+1,polyval(polyfit(actual(keep),predictors(keep,2)',1),unique(actual(keep))))
title([num2str(i) ', mu = ' num2str(mseZ_cells((idx(i)),2)) ', R=' num2str(corrs(idx(i),2)) ', p=' num2str(pval(idx(i),2))])
ylabel(names{2})
xlabel('LS # of spks')

subplot(2,2,2)
% scatter(actual(keep),predictors(keep,3),'.k')
 boxplot(predictors(keep,3),actual(keep),'color','k')
hold on
plot(unique(actual(keep))+1,polyval(polyfit(actual(keep),predictors(keep,3)',1),unique(actual(keep))))
[c p] = corr(actual(keep)',predictors(keep,3));
title(['mu = ' num2str(mseZ_cells((idx(i)),3)) ', R=' num2str(c) ', p=' num2str(p)])
ylabel(names{3})
xlabel('LS # of spks')

subplot(2,2,3)
% scatter(actual(keep),predictors(keep,9),'.k')
 boxplot(predictors(keep,9),actual(keep),'color','k')
hold on
% keep2 = find(~isnan((predictors(:,9)))); 
plot(unique(actual(keep2))+1,polyval(polyfit(actual(keep2),predictors(keep2,9)',1),unique(actual(keep2))))
title(['mu = ' num2str(mseZ_cells((idx(i)),9)) ', R=' num2str(corrs(idx(i),9)) ', p=' num2str(pval(idx(i),9))])
ylabel(names{9})
xlabel('LS # of spks')

subplot(2,2,4)
% scatter(actual(keep),predictors(keep,8),'.k')
 boxplot(predictors(keep,8),actual(keep),'color','k')
hold on
% keep2 = find(~isnan((predictors(:,8)))); 
plot(unique(actual(keep2))+1,polyval(polyfit(actual(keep2),predictors(keep2,8)',1),unique(actual(keep2))))
title(['mu = ' num2str(mseZ_cells((idx(i)),8)) ', R=' num2str(corrs(idx(i),8)) ', p=' num2str(pval(idx(i),8))])
ylabel(names{8})
xlabel('LS # of spks')

pause
end
end