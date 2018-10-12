% 
clear all
d = dir('*201*');
for i=1:length(d)
cd(d(i).name)
if exist([d(i).name '.content_GLM.mat']) & exist([d(i).name '.bayesianResults_popBursts.mat']) % or ripple bayesianResults
dat = load([d(i).name '.content_GLM.mat']);
content(i) = dat.content;
dat = load([d(i).name '.bayesianResults_popBursts.mat']);
if isfield(dat,'integral_hpc')
bay{i} = dat.integral_hpc;
else
bay{i} = nan(2,length(content(i).nSpikes{1}));
end
end
cd ~/datasets/ripples_LS/
end

count = 1;

for rec = 1:length(content)
    if ~isempty(content(rec).nSpikes)
        idx_hpc = find(strcmp(content(rec).region{1},'hpc'));
        if isempty(idx_hpc)
            idx_hpc = find(strcmp(content(rec).region{1},'ca3'));
        end
        
predictors = [nanmean(content(rec).nSpikes{1}(idx_hpc,:))' ...
              content(rec).hpc_popBurst{1}...
              [(nanmean(content(rec).PF{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0)))./nanmean(content(rec).nSpikes{1}(idx_hpc,:))]'... % % of spikes from PF
              [[nanmean(content(rec).spatialContent{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0))] ./nanmean(content(rec).nSpikes{1}(idx_hpc,:))]' ... 
              [[nanmean(content(rec).rewardContent{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0))] ./ nanmean(content(rec).nSpikes{1}(idx_hpc,:))]'...
              nanmean(content(rec).cellLoc_wav{1}(idx_hpc,:).*content(rec).nSpikes{1}(idx_hpc,:))'...
              content(rec).condition{1}...
              content(rec).location{1}...
              content(rec).SleepState{1}...
              content(rec).duration{1}'...
              nanmax(bay{rec})'];%...
%               bay{rec}'];
          
actual = content(rec).ls_popBurst{1};
% actual = content(rec).ls_power{1}';

for p = 1:size(predictors,2)
    if ~isnan(nansum(predictors(:,p)))
        [beta dev(p) ] = glmfit(predictors(:,p),actual,'normal');
        yfit = glmval(beta,predictors(:,p),'identity');
        mse(p) = nanmean((yfit-actual').^2);

        for iter = 1:100
            pred_shuf = bz_shuffleCircular(predictors(:,p)')';
%             pred_shuf = predictors(randperm(size(predictors,1)),p);
            [beta dev_shuf(p,iter)] = glmfit(pred_shuf,actual,'normal');
            yfit = glmval(beta,pred_shuf,'identity');
            mse_shuf(p,iter) = nanmean((yfit-actual').^2);
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
boundedline(1:size(devZ,2),nanmean(devZ),sem(devZ))



%%
keep = find(~isnan(predictors(:,end)));
%%
if ~isempty(keep)


idx = find(strcmp(content(rec).region{1},'ls'));
hpc = find(strcmp(content(rec).region{1},'hpc'));
if ~isempty(hpc) 
    hpc = find(strcmp(content(rec).region{1},'ca3'));
end

for spk = 1:length(idx)
   actual = content(rec).nSpikes{1}(idx(spk),keep); 
   if mean(actual) > 0
       for p = 1:size(predictors,2)
        if ~isnan(nansum(predictors(:,p)))
            [beta dev(p) ] = glmfit(predictors(keep,p),actual,'normal');
            yfit = glmval(beta,predictors(keep,p),'identity');
            mse(p) = nanmean((yfit-actual').^2);

            for iter = 1:100
                pred_shuf = bz_shuffleCircular(predictors(keep,p)')';
    %             pred_shuf = predictors(randperm(size(predictors,1)),p);
                [beta dev_shuf(p,iter)] = glmfit(pred_shuf,actual,'normal');
                yfit = glmval(beta,pred_shuf,'identity');
                mse_shuf(p,iter) = nanmean((yfit-actual').^2);
            end
        else
            dev(p) = nan;
            dev_shuf(p,1:100) = nan;
        end
       end
       
       meanCount(count) = mean(actual);
       mseZ_cells(count,:) = (mse-nanmean(mse_shuf,2)')./nanstd(mse_shuf,[],2)';
       mse_cells(count,:) = mse;
       mse_shuf_cells(count,:,:) = mse_shuf;
       recording(count) = rec;
       count = 1+count;
   end
end


% subplot(2,2,4)
% cla
% boundedline(1:size(mseZ_cells,2),nanmean(mseZ_cells),sem(mseZ_cells))

    end


pause(.001)
    end
end
%% CHECK THAT ls BURSTS ARE LARGER AFTER BEHAV????
%  for i=1:3
% idx = find(content(rec).condition{1}==i);
% plot(i,nanmean(content(rec).ls_popBurst{1}(idx)),'.r')
% hold on
% end