% 
clear all
d = dir('*201*');
% for i=1:length(d)
% cd(d(i).name)
% if exist([d(i).name '.content_GLM.mat']) & ...
%         exist([d(i).name '.bayesianResults_popBursts.mat']) & ...
%         exist([d(i).name '.rankOrder_popBursts.mat'])% & ...
% %         exist([d(i).name '_seqNMF.popBursts.hpc.mat']) 
%     
% dat = load([d(i).name '.content_GLM.mat']);
% content(i) = dat.content;
% 
% dat = load([d(i).name '.bayesianResults_popBursts.mat']);
% if isfield(dat,'integral_hpc')
% bay{i} = dat.integral_hpc;
% else
% bay{i} = nan(2,length(content(i).nSpikes{1}));
% end
% 
% dat = load([d(i).name '.rankOrder_popBursts.mat']);
% rankOrder{i} = dat.rankOrder;
% 
% %% seqNMF data here
% % seq = load([d(i).name '_seqNMF.popBursts.hpc.mat']);
% % for e=1:length(dat.integral_hpc)
% %     ts= e*101+49;
% %     if ts+25<=length(seq.H_hpc)
% %         [vec] = max(seq.H_hpc(:,ts-25:ts+25)');
% %     else
% %         vec = nan; 
% %     end
% %     seqs{i}(e,:) = vec;
% % end
% 
% end
% cd ~/datasets/ripples_LS/
% end


for i=1:length(d)
cd(d(i).name)
if exist([d(i).name '.content_GLM_ripple.mat']) & ...
        exist([d(i).name '.bayesianResults_ripples.mat'])& ...
        exist([d(i).name '.rankOrder_ripples.mat'])% & ...
%         exist([d(i).name '_seqNMF.ripples.hpc.mat']) 
    
dat = load([d(i).name '.content_GLM_ripple.mat']);
content(i) = dat.content;

dat = load([d(i).name '.bayesianResults_ripples.mat']);
if isfield(dat,'integral_hpc')
bay{i} = dat.integral_hpc;
else
bay{i} = nan(2,length(content(i).nSpikes{1}));
end

dat = load([d(i).name '.rankOrder_ripples.mat']);
rankOrder{i} = dat.rankOrder;

%% seqNMF data here
% seq = load([d(i).name '_seqNMF.ripples.hpc.mat']);
% for e=1:length(dat.integral_hpc)
%     ts= e*101+49;
%     if ts+25<=length(seq.H_hpc)
%         [vec] = max(seq.H_hpc(:,ts-25:ts+25)');
%     else
%         vec = nan; 
%     end
%     seqs{i}(e,:) = vec;
% end

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
        
        names = {'# spikes','PF spk %',...
                'spatial content',...
                'reward content',...
                'ripple depth',...
                'pre/behav/post condition',...
                'location',...
                'Sleep state',....
                'duration',...
                'rs1 (radon int)',...
                'rs2 (rank ord)'};
            
            
predictors = [nanmean(content(rec).nSpikes{1}(idx_hpc,:))' ...
              [(nanmean(content(rec).PF{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0)))./nanmean(content(rec).nSpikes{1}(idx_hpc,:))]'... % % of spikes from PF
              [[nanmean(content(rec).spatialContent{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0))] ./nanmean(content(rec).nSpikes{1}(idx_hpc,:))]' ... 
              [[nanmean(content(rec).rewardContent{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0))] ./ nanmean(content(rec).nSpikes{1}(idx_hpc,:))]'...
              nanmean(content(rec).cellLoc_wav{1}(idx_hpc,:).*content(rec).nSpikes{1}(idx_hpc,:))'...
              content(rec).condition{1}...
              content(rec).location{1}...
              content(rec).SleepState{1}...
              content(rec).duration{1}'...
              nanmax(bay{rec})'];%...
%               nanmean(rankOrder{rec})'];%...
%               bay{rec}'];
%               content(rec).hpc_popBurst{1}...
          
actual = content(rec).ls_popBurst{1};
% actual = content(rec).ls_power{1}';

for p = 1:size(predictors,2)
    if ~isnan(nansum(predictors(:,p)))
        [beta dev(p) ] = glmfit([ predictors(:,p)],actual','normal');
        yfit = glmval(beta,[predictors(:,p)],'identity');
        mse(p) = nanmean((yfit-actual').^2);

        for iter = 1:100
%             pred_shuf = [predictors(:,1:p-1) bz_shuffleCircular(predictors(:,p)')' predictors(:,p+1:end) ];
            pred_shuf = [bz_shuffleCircular(predictors(:,p)')'  ];
            
%             pred_shuf = predictors(randperm(size(predictors,1)),p);
            [beta dev_shuf(p,iter)] = glmfit([ pred_shuf],actual','normal');
            yfit = glmval(beta,[ pred_shuf],'identity');
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

%%
idx = find(strcmp(content(rec).region{1},'ls'));
hpc = find(strcmp(content(rec).region{1},'hpc'));
if ~isempty(hpc) 
    hpc = find(strcmp(content(rec).region{1},'ca3'));
end

% keep = find(~isnan(sum(predictors,2))); 

for spk = 1:length(idx)
   
   if mean(actual) > 0
       for p = 1:size(predictors,2)
           keep = find(~isnan(sum(predictors(:,p)))); 
           actual = content(rec).nSpikes{1}(idx(spk),keep); 
        if ~isnan(nansum(predictors(:,p))) & length(keep) > 100
            [beta dev(p) ] = glmfit(predictors(keep,p),actual,'normal');
            yfit = glmval(beta,predictors(keep,p),'identity');
            mse(p) = nanmean((yfit-actual').^2);

            
            [corrs(count,p)] = corr(predictors(keep,p),actual','rows','complete');
       
            for iter = 1:100
%                  pred_shuf = [predictors(keep,1:p-1) bz_shuffleCircular(predictors(keep,p)')' predictors(keep,p+1:end) ];
                pred_shuf = bz_shuffleCircular(predictors(keep,p)')';
    %             pred_shuf = predictors(randperm(size(predictors,1)),p);
                [beta dev_shuf(p,iter)] = glmfit(pred_shuf,actual,'normal');
                yfit = glmval(beta,pred_shuf,'identity');
                mse_shuf(p,iter) = nanmean((yfit-actual').^2);
                
                [corrs_shuff(count,p,iter)] = corr(pred_shuf,actual','rows','complete');
            end
                    
        [ft r1]= fit(predictors(keep,p),actual','poly1');
        [ft r2]= fit(predictors(keep,p),actual','poly2');
        adj_R(p,:) = [r1.adjrsquare r2.adjrsquare];
        
        else
            corrs(count,p) = nan;
            dev(p) = nan;
            dev_shuf(p,1:100) = nan;
            adj_R(p,1:2) = nan;
        end
       end
       
       meanCount(count) = nanmean(actual);
       mseZ_cells(count,:) = (mse-nanmean(mse_shuf,2)')./nanstd(mse_shuf,[],2)';
       adjRSqaure(count,:,:) = adj_R;
%        mseZ_cells(count,:) = (mse_shuf(:,:,1)-nanmean(mse_shuf,2)')./nanstd(mse_shuf,[],2)';
       mse_cells(count,:) = mse;
       mse_shuf_cells(count,:,:) = mse_shuf;
       recording(count) = rec;
       count = 1+count;
   end
end


subplot(2,2,4)
cla
boundedline(1:size(mseZ_cells,2),nanmean(mseZ_cells,1),sem(mseZ_cells))


pause(.001)
    end
end
%% CHECK THAT ls BURSTS ARE LARGER AFTER BEHAV????
%  for i=1:3
% idx = find(content(rec).condition{1}==i);
% plot(i,nanmean(content(rec).ls_popBurst{1}(idx)),'.r')
% hold on
% end

idx = find(meanCount>.1);
[a b] = sort(nanmean(mseZ_cells(idx,:)));

figure

for i=1:size(mseZ_cells,2)
subplot(size(mseZ_cells,2),1,i)
raincloud_plot('X',corrs(idx,i),'density_type', 'ks','bandwidth',.051,'color',[0 0 0]);
raincloud_plot('X',corrs_shuff(idx,i,1),'density_type', 'ks','bandwidth',.051,'color',[1 0 0]);


% raincloud_plot('X',mseZ_cells(idx,i),'density_type', 'ks','bandwidth',.2,'color',[0 0 0]);
% raincloud_plot('X',mseZ_cells_shuf(idx,i,1),'density_type', 'ks','bandwidth',.2,'color',[1 0 0]);
% axis([-20 2 -1 1])
title(names{i})
end


