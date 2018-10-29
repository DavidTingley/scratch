idx = find(pval(:,1)<.05 & corrs(:,1) > .1);

for i=1:length(idx)
rec = recording(idx(i));
idx_hpc = find(strcmp(content(rec).region{1},'hpc'));
if isempty(idx_hpc)
idx_hpc = find(strcmp(content(rec).region{1},'ca3'));
end
predictors = [nansum(content(rec).nSpikes{1}(idx_hpc,:))' ...
[(nanmean(content(rec).PF{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0)))./nanmean(content(rec).nSpikes{1}(idx_hpc,:)~=0)]'... % % of spikes from PF
[[nanmean(content(rec).rewardContent{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)))] ./ nanmean(content(rec).nSpikes{1}(idx_hpc,:))]'...
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
scatter(actual(keep),predictors(keep,1),'.k')
hold on
plot(unique(actual(keep)),polyval(polyfit(actual(keep),predictors(keep,1)',1),unique(actual(keep))))
title(['R=' num2str(corrs(idx(i),1)) ', p=' num2str(pval(idx(i),1))])
ylabel(names{1})
subplot(2,2,2)
scatter(actual(keep),predictors(keep,7),'.k')
hold on
plot(unique(actual(keep)),polyval(polyfit(actual(keep),predictors(keep,7)',1),unique(actual(keep))))
title(['R=' num2str(corrs(idx(i),7)) ', p=' num2str(pval(idx(i),7))])
ylabel(names{7})

subplot(2,2,3)
scatter(actual(keep),predictors(keep,9),'.k')
hold on
% keep2 = find(~isnan((predictors(:,9)))); 
plot(unique(actual(keep2)),polyval(polyfit(actual(keep2),predictors(keep2,9)',1),unique(actual(keep2))))
title(['R=' num2str(corrs(idx(i),9)) ', p=' num2str(pval(idx(i),9))])
ylabel(names{9})

subplot(2,2,4)
scatter(actual(keep),predictors(keep,8),'.k')
hold on
% keep2 = find(~isnan((predictors(:,8)))); 
plot(unique(actual(keep2)),polyval(polyfit(actual(keep2),predictors(keep2,8)',1),unique(actual(keep2))))
title(['R=' num2str(corrs(idx(i),8)) ', p=' num2str(pval(idx(i),8))])
ylabel(names{8})
pause
end
end