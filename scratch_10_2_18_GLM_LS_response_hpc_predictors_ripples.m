% clear all
d = dir('*201*');
warning off
cd  ~/datasets/ripples_LS/

for i=1:length(d)
cd(d(i).name)
if exist([d(i).name '.content_GLM_ripple_25ms.mat'])% & ...
%         exist([d(i).name '.bayesianResults_popBursts.mat']) & ...
%         exist([d(i).name '_seqNMF.popBursts.hpc.mat']) 
%         exist([d(i).name '.rankOrder_popBursts.mat'])% & ...
        
    
dat = load([d(i).name '.content_GLM_ripple_25ms.mat'],'content');
content(i) = dat.content;
% 
if  exist([d(i).name '.replayComparisonData_oneField.mat']) 
    dat = load([d(i).name '.replayComparisonData_oneField.mat'],'*slope*', 'integral_hpc*','corr*','outR*','rankOrder');
    if isfield(dat,'integral_hpc') & isfield(dat,'outR')
    radonInt{i} = dat.integral_hpc;
    radonSlope{i} = dat.slope_hpc;
    weightedCorr{i} = dat.outR;
    replay{i}.rZscore = (abs(dat.outR) - nanmean(abs(dat.outR_shuf),3)) ./ nanstd(abs(dat.outR_shuf),[],3);
    rankOrder{i} = absmax(dat.rankOrder);
    replay{i}.rankOrder = dat.rankOrder;
    
     for j = 1:size(dat.corrs_hpc,1)
        for k=1:size(dat.corrs_hpc,2)
%             [a replay{i}.rankOrder_bayes_pvals(j,k)] = ttest2(dat.corrs_hpc(j,k),dat.corrs_hpc_shuf(j,k,:));
            [a replay{i}.weightedCorr_pvals(j,k)] = ttest2(dat.outR(j,k),dat.outR_shuf(j,k,:));
        end
    end
    replay{i}.radonIntegral = dat.integral_hpc;
    replay{i}.radonSlope = dat.slope_hpc;
    replay{i}.radonIntegral_z = radonInt{i};
    replay{i}.linearWeighted = weightedCorr{i};
    
    else
    radonInt{i} = nan(2,size(content(i).nSpikes{1},2));
    weightedCorr{i} = nan(2,size(content(i).nSpikes{1},2));
    end
else
    radonInt{i} = nan(2,size(content(i).nSpikes{1},2));
    weightedCorr{i} = nan(2,size(content(i).nSpikes{1},2));
    rankOrder{i} = nan(1,size(content(i).nSpikes{1},2));
end

%% rank order correlations here
if isfield(dat,'corrs_hpc')
    rankOrder_bayes{i} = absmax(dat.corrs_hpc); %dat.rankOrder;
    for j = 1:size(dat.corrs_hpc,1)
        for k=1:size(dat.corrs_hpc,2)
            [a replay{i}.rankOrder_bayes_pvals(j,k)] = ttest2(dat.corrs_hpc(j,k),dat.corrs_hpc_shuf(j,k,:));
%             [a replay{i}.rankOrder_pvals(j,k)] = ttest2(dat.outR(j,k),dat.outR_shuf(j,k,:));
        end
    end
    
    replay{i}.rankOrder_bayes =  dat.corrs_hpc;
else
    rankOrder_bayes{i} = nan(1,size(content(i).nSpikes{1},2));
end






%% seqNMF data here
if exist([d(i).name '_seqNMF.ripples.hpc.mat']) 
seq = load([d(i).name '_seqNMF.ripples.hpc.mat']);
if ~isempty(content(i))
    for e=1:size(content(i).nSpikes{1},2)
        ts= e*101+49;
        if ts+25<=length(seq.H_hpc)
            [vec] = max(seq.H_hpc(:,ts-25:ts+25)');
        else
            vec = nan; 
        end
        seqs{i}(e,:) = vec;
    end
else
    seqs{i} = nan(size(content(i).nSpikes{1},2),2);
end
else
     seqs{i} = nan(size(content(i).nSpikes{1},2),2);
end
end


%% only coupled events, reviewer #1              -_____-
ripples = bz_LoadEvents(pwd,'CA1Ripples');
hfo = bz_LoadEvents(pwd,'LSRipples');
if ~isempty(ripples) & ~isempty(hfo)
    for t = 1:length(ripples.peaks)
        coupled{i}(t) = 0;
        [a b] = min(abs(ripples.peaks(t)-hfo.peaks));
        if ripples.peaks(t) - hfo.peaks(b) < 0 & a > .05
           coupled{i}(t) = 1; 
        end
    end
end
% cd .. 
% cd D:\datasets\ripples_LS
cd ~/datasets/ripples_LS/
end


% for i=1:length(d)
% cd(d(i).name)
% if exist([d(i).name '.content_GLM_ripple_25ms.mat']) & ...
%         exist([d(i).name '.bayesianResults_ripples.mat'])% & ...
% %         exist([d(i).name '.rankOrder_ripples.mat'])% & ...
% %         exist([d(i).name '_seqNMF.ripples.hpc.mat']) 
%     
% dat = load([d(i).name '.content_GLM_ripple_25ms.mat'],'content');
% content(i) = dat.content;
% 
% dat = load([d(i).name '.bayesianResults_ripples.mat'],'integral_hpc*','corr*');
% if isfield(dat,'integral_hpc')
% radonInt{i} = absmax((dat.integral_hpc - nanmean(dat.integral_hpc_shuf,3)) ./ nanstd(dat.integral_hpc_shuf,[],3));
% else
% radonInt{i} = nan(2,length(content(i).nSpikes{1}));
% end
% % 
% % dat = load([d(i).name '.rankOrder_ripples.mat'],'rankOrder');
% rankOrder{i} = absmax(dat.corrs_hpc); %dat.rankOrder;
% 
% %% seqNMF data here
% % seq = load([d(i).name '_seqNMF.ripples.hpc.mat']);
% % for e=1:length(dat.rankOrder)
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



count = 1;

for rec = 1:length(content)
    if ~isempty(content(rec).nSpikes)
        idx_hpc = find(strcmp(content(rec).region{1},'hpc'));
        if isempty(idx_hpc)
            idx_hpc = find(strcmp(content(rec).region{1},'ca3'));
        end
        
        names = {'# spikes',...%                  '% PF spks',...%                 'reward content',...
                'inter-ripple interval',...
                'ripple depth',...
                'pre/behav/post condition',...
                'location',...
                'Sleep state',....
                'duration',...
                'rs1 (radon int)',...
                'rs2 (rank ord,bayes)',...
                'rs3 (rank ord, raw)',...
                'rs4 (weighted Corr)',...
                'seqNMF',...
                'hpc cellID vec',...
                'seqNMF_all'};
            
for i = 1:size(seqs{rec},1)
    [seq_fits(i) seq_ids(i)] = max(seqs{rec}(i,:));
end

i
predictors = [nanmean(content(rec).nSpikes{1}(idx_hpc,:))' ... %               [(nanmean(content(rec).PF{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)~=0)))./nanmean(content(rec).nSpikes{1}(idx_hpc,:)~=0)]'... % % of spikes from PF%               [[nanmean(content(rec).rewardContent{1}(idx_hpc,:) .* (content(rec).nSpikes{1}(idx_hpc,:)))] ./ nanmean(content(rec).nSpikes{1}(idx_hpc,:))]'...
              content(rec).interRipInterval{1}',...
              nanmean(content(rec).cellLoc_wav{1}(idx_hpc,:).*content(rec).nSpikes{1}(idx_hpc,:))'...
              content(rec).condition{1}...
              content(rec).location{1}...
              content(rec).SleepState{1}...
              content(rec).duration{1}'...
              max(radonInt{rec})'...
              rankOrder_bayes{rec}'...
              rankOrder{rec}'...
              absmax(weightedCorr{rec})'...
              max(seqs{rec}')'];%...
          
% actual = content(rec).ls_popBurst{1};
% actual = content(rec).ls_power{1}';

% for p = 1:size(predictors,2)
%     if ~isnan(nansum(predictors(:,p)))
%         [beta dev(p) ] = glmfit([ predictors(:,p)],actual','normal');
%         yfit = glmval(beta,[predictors(:,p)],'identity');
%         mse(p) = nanmean((yfit-actual').^2);
% 
%         for iter = 1:100
% %             pred_shuf = [predictors(:,1:p-1) bz_shuffleCircular(predictors(:,p)')' predictors(:,p+1:end) ];
%             pred_shuf = [bz_shuffleCircular(predictors(:,p)')'  ];
%             
% %             pred_shuf = predictors(randperm(size(predictors,1)),p);
%             [beta dev_shuf(p,iter)] = glmfit([ pred_shuf],actual','normal');
%             yfit = glmval(beta,[ pred_shuf],'identity');
%             mse_shuf(p,iter) = nanmean((yfit-actual').^2);
%         end
% 
%     else
%         dev(p) = nan;
%         dev_shuf(p,1:100) = nan;
%     end
% end
% 
% 
% subplot(2,2,1)
% cla
% boundedline(1:size(predictors,2),nanmean(dev_shuf,2),std(dev_shuf,[],2))
% plot(dev,'r')
% subplot(2,2,2)
% devZ(rec,:) = (dev-nanmean(dev_shuf,2)')./std(dev_shuf,[],2)';
% devs(rec,:) = dev; clear dev dev_shuf
% imagesc(devZ)
% subplot(2,2,3)
% cla
% boundedline(1:size(devZ,2),nanmean(devZ),sem(devZ))



%%

%%
idx = find(strcmp(content(rec).region{1},'ls'));
hpc = find(strcmp(content(rec).region{1},'hpc'));
if isempty(hpc) 
    hpc = find(strcmp(content(rec).region{1},'ca3'));
end

% keep = find(~isnan(sum(predictors,2))); 
s = seqs{rec};

for spk = 1:length(idx)
%       if mean(actual) > 0
       %% now for most predictors
       counts = sum(double(content(rec).nSpikes{1}(hpc,:)))'; 
%        counts = sum(double(content(rec).nSpikes{1}(hpc,:)))'; 

       for p = 1:size(predictors,2)
           keep = find(~isnan((predictors(:,p)))); 
           %% for excluding 'non-coupled' events.
           keep = intersect(keep,find(coupled{rec}));
           
           actual = content(rec).nSpikes{1}(idx(spk),keep); 
        if ~isnan(nansum(predictors(:,p))) & length(keep) > 10
            [beta dev(p) ] = glmfit([predictors(keep,p)],actual,'normal');
            yfit = glmval(beta,[predictors(keep,p)],'identity');
%             [beta dev(p) ] = glmfit(predictors(keep,:),actual,'normal');
%             yfit = glmval(beta,predictors(keep,:),'identity');
            mse(p) = nanmean((yfit-actual').^2);            
            [corrs(count,p) pval(count,p)] = corr(predictors(keep,p),actual','rows','complete','type','spearman');
       
            for iter = 1:100
%                  pred_shuf = [predictors(keep,1:p-1) bz_shuffleCircular(predictors(keep,p)')' predictors(keep,p+1:end) ];
                pred_shuf = bz_shuffleCircular(predictors(keep,p)')';
    %             pred_shuf = predictors(randperm(size(predictors,1)),p);
                [beta dev_shuf(p,iter)] = glmfit([pred_shuf],actual,'normal');
                yfit = glmval(beta,[pred_shuf],'identity');
                mse_shuf(p,iter) = nanmean((yfit-actual').^2);
                
                

                
%                 [corrs_shuff(count,p,iter)] = corr(pred_shuf(:,p),actual','rows','complete');
                [corrs_shuff(count,p,iter) pval_shuff(count,p,iter)] = corr(pred_shuf,actual','rows','complete','type','spearman');
            end
                    
        [ft r1]= fit(predictors(keep,p),actual','poly1');
        [ft r2]= fit(predictors(keep,p),actual','poly2');
        adj_R(p,:) = [r1.adjrsquare r2.adjrsquare];
        
        else
            mse(p) = nan;
            mse_shuf(p,1:100) = nan;
            corrs(count,p) = nan;
            pval(count,p) = nan;
            pval_shuff(count,p,1:100) = nan;
            dev(p) = nan;
            dev_shuf(p,1:100) = nan;
            adj_R(p,1:2) = nan;
        end
       end
       
       %% now for cell ID vec prediction
       actual = content(rec).nSpikes{1}(idx(spk),:); 
       pred = double(content(rec).nSpikes{1}(hpc,:)>0)'; 
       counts = sum(double(content(rec).nSpikes{1}(hpc,:)))'; 
%        counts = sum(double(content(rec).nSpikes{1}(hpc,:)))'; 
       [beta dev_cellID stats{count}] = glmfit([pred mean(pred')' counts],actual','normal');
       yfit = glmval(beta,[pred mean(pred')' counts],'identity');
       mse_cellIDVec = nanmean((yfit-actual').^2);
       clear pred_shuf
%        [beta dev_pop] = glmfit([counts],actual','normal');
       for cell = 1:size(pred,2)
           [beta dev_cell(cell)] = glmfit([pred(:,cell) counts],actual','normal');
       end
       
       for iter = 1:100
%            for p = 1:size(pred,1)
%            pred_shuf(p,:) = bz_shuffleCellID(pred(p,:)');
%            end
           pred_shuf = bz_shuffleCircular(pred);
           
           [beta dev_cellID_shuf(iter)] = glmfit([pred_shuf mean(pred_shuf')' counts],actual','normal');
           yfit = glmval(beta,[pred_shuf mean(pred_shuf')' counts],'identity');
           mse_cellIDVec_shuf(iter) = nanmean((yfit-actual').^2);
           
           for cell = 1:size(pred,2)
                [beta dev_cell_shuf(iter,cell)] = glmfit([pred_shuf(:,cell) counts],actual','normal');
           end
       
           
       end
       mse_individual_cells{count} = (dev_cell-mean(dev_cell_shuf))./nanstd(dev_cell_shuf); 
       for iter = 1:100
            mse_individual_cells_shuf{count}(iter,:) = (mean(dev_cell_shuf)-mean(mean(dev_cell_shuf)))./nanstd(mean(dev_cell_shuf)); 
       end
       clear dev_cell* dev_pop
       mseZ_cells(count,13) = (mse_cellIDVec - mean(mse_cellIDVec_shuf,2)) ./ nanstd(mse_cellIDVec_shuf,[],2);
       mseZ_cells_shuf(count,13) = (mse_cellIDVec_shuf(1) - mean(mse_cellIDVec_shuf(2:end),2)) ./ nanstd(mse_cellIDVec_shuf(2:end),[],2);
%        pvals_cellID(count,:) = stats.p(3:end);
       
       %% now for seqNMF prediction
       actual = content(rec).nSpikes{1}(idx(spk),:); 
       [beta dev_seqNMF] = glmfit([predictors(:,1) (s)],actual','normal');
       yfit = glmval(beta,[predictors(:,1) s],'identity');
       mse_seqNMF = nanmean((yfit-actual').^2);
       
       for iter = 1:100
             % now for seqNMF
            s_shuf = bz_shuffleCircular(s);
            [beta dev_seqNMF_shuf(iter)] = glmfit([predictors(:,1) s_shuf],actual','normal');  
            yfit = glmval(beta,[predictors(:,1) s_shuf],'identity');
            mse_seqNMF_shuf(iter) = nanmean((yfit-actual').^2);
       end
       
       mseZ_cells(count,14) = (mse_seqNMF - nanmean(mse_seqNMF_shuf)) ./ nanstd(mse_seqNMF_shuf);
       mseZ_cells_shuf(count,14) = (mse_seqNMF_shuf(1) - nanmean(mse_seqNMF_shuf(2:end))) ./ nanstd(mse_seqNMF_shuf(2:end));
       
       meanCount(count) = nanmean(actual);
       mseZ_cells(count,1:12) = (mse-nanmean(mse_shuf,2)')./nanstd(mse_shuf,[],2)';
       adjRSqaure(count,:,:) = adj_R;
       mseZ_cells_shuf(count,1:12) = (squeeze(mse_shuf(:,1))-nanmean(mse_shuf,2))./nanstd(mse_shuf,[],2);
       mse_cells(count,:) = mse;
       mse_shuf_cells(count,:,:) = mse_shuf;
       seqNMF_mse(count) = mse_seqNMF;
       seqNMF_mse_shuf(count,:) = mse_seqNMF_shuf;
       
       recording(count) = rec;
       cellID(count) = idx(spk);
       region{count} = content(rec).region{1}{idx(spk)};
       count = 1+count;
%    end
end


subplot(2,2,4)
cla
boundedline(1:size(mseZ_cells,2),nanmean(mseZ_cells,1),sem(mseZ_cells))


pause(.001)
    end
end
%% CHECK THAT ls BURSTS ARE LARGER AFTER BEHAV????
% for i=1:3
% idx = find(content(rec).condition{1}==i);
% plot(i,nanmean(content(rec).ls_popBurst{1}(idx)),'.r')
% hold on
% end

mseZ_cells(mseZ_cells<-1000) = nan;
idx = find(meanCount>.05);
nCellCount = sum(~isnan(mseZ_cells));
[a b] = sort(nanmean(mseZ_cells(idx,:)),'descend');
% [a b] = sort(nanmean(pval(idx,:)<.05),'descend');

figure

for i=1:size(mseZ_cells,2)
subplot(size(mseZ_cells,2),1,i)
% raincloud_plot('X',corrs(idx,b(i)),'density_type', 'ks','bandwidth',.025,'color',[0 0 0]);
% raincloud_plot('X',corrs_shuff(idx,b(i),1),'density_type', 'ks','bandwidth',.025,'color',[1 0 0]);

raincloud_plot('X',mseZ_cells(:,b(i)),'density_type', 'ks','bandwidth',.2,'color',[0 0 0]);
raincloud_plot('X',mseZ_cells_shuf(:,b(i),1),'density_type', 'ks','bandwidth',.2,'color',[1 0 0]);
% axis([-110 2 -1 3])
xlim([-150 3])
title([names{b(i)} ', ' num2str(sum(mseZ_cells(:,b(i))<-2)./sum(~isnan(mseZ_cells(:,b(i))))) '% P < .05, mean=' num2str(nanmean(mseZ_cells(:,b(i)))) ', N=' num2str(nCellCount(b(i)))])

end


figure
for i=1:size(mseZ_cells,2)
m{i,1} =mseZ_cells(:,i);
m{i,2} =mseZ_cells_shuf(:,i);
end
clz{1} = repmat([0 0 0],14,1);
clz{2} = repmat([1 0 0],14,1);
raincloud_lineplot_2(m,clz,1,1)


save('C:\Users\SB13FLLT001\Dropbox\Rip_features_seqNMF_LS_data_ripples_25ms.mat','-v7.3')

