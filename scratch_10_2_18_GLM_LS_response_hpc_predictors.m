
% rec = 4;
for rec = 1:97
    if ~isempty(content.nSpikes{rec})
predictors = [nanmean(content.nSpikes{rec})' ...
              nanmean(content.PF{rec})'...
              nanmean(content.spatialContent{rec})' ... % check that this isn't just FR
              nanmean(content.rewardContent{rec})'...
              nanmean(content.cellLoc_wav{rec}.*content.nSpikes{rec})'...
              content.condition{rec}...
              content.location{rec}...
              content.SleepState{rec}...
              double(content.hpc_popBurst{rec})];
          
% actual = content.ls_popBurst{rec};
actual = content.ls_power{rec}';

for p = 1:size(predictors,2)
    if ~isnan(sum(predictors(:,p)))
        [beta dev(p) ] = glmfit(predictors(:,p),actual,'normal');
        yfit = glmval(beta,predictors(:,p),'identity');
        sse(p) = sum((yfit-actual').^2);

        for iter = 1:100
            pred_shuf = bz_shuffleCircular(predictors(:,p)')';
    %         pred_shuf = predictors(randperm(size(predictors,1)),p);
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
imagesc(devZ)
pause(.001)
    end
end
%% CHECK THAT ls BURSTS ARE LARGER AFTER BEHAV????
 for i=1:3
idx = find(content.condition{rec}==i);
plot(i,nanmean(content.ls_popBurst{rec}(idx)),'.r')
hold on
end