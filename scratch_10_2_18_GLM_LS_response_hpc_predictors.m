
% rec = 4;
for rec = 1:97
    if ~isempty(content(rec).nSpikes)
predictors = [nanmean(content(rec).nSpikes{1})' ...
              nanmean(content(rec).PF{1})'...
              nanmean(content(rec).spatialContent{1})' ... % check that this isn't just FR
              nanmean(content(rec).rewardContent{1})'...
              nanmean(content(rec).cellLoc_wav{1}.*content(rec).nSpikes{1})'...
              content(rec).condition{1}...
              content(rec).location{1}...
              content(rec).SleepState{1}...
              double(content(rec).hpc_popBurst{1})];
          
actual = content(rec).ls_popBurst{1};
% actual = content(rec).ls_power{1}';

for p = 1:size(predictors,2)
    if ~isnan(sum(predictors(:,p)))
        [beta dev(p) ] = glmfit(predictors(:,p),actual,'normal');
        yfit = glmval(beta,predictors(:,p),'identity');
        sse(p) = sum((yfit-actual').^2);

        for iter = 1:100
            pred_shuf = bz_shuffleCircular(predictors(:,p)')';
%             pred_shuf = predictors(randperm(size(predictors,1)),p);
            [beta dev_shuf(p,iter)] = glmfit(pred_shuf,actual,'normal');
            yfit = glmval(beta,pred_shuf,'identity');
            sse_shuf(p,iter) = sum((yfit-actual').^2);
        end
    else
        dev(p) = nan;
        dev_shuf(p,1:100) = nan;
    end
end


subplot(2,2,1)
cla
boundedline(1:size(predictors,2),nanmean(dev_shuf,2),std(dev_shuf,[],2))
plot(dev,'r')
subplot(2,2,2)
devZ(rec,:) = (dev-nanmean(dev_shuf,2)')./std(dev_shuf,[],2)';
devs(rec,:) = dev; clear dev dev_shuf
imagesc(devZ)
subplot(2,2,3)
cla
boundedline(1:size(devZ,2),nanmedian(devZ),sem(devZ))
pause(.001)
    end
end
%% CHECK THAT ls BURSTS ARE LARGER AFTER BEHAV????
 for i=1:3
idx = find(content(rec).condition{1}==i);
plot(i,nanmean(content(rec).ls_popBurst{1}(idx)),'.r')
hold on
end