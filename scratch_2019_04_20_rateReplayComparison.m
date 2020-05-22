binSize = .01;
for c = 5:100
    for iter =1:100
    iter
    for i=1:10
        r = randperm(100);
        data = data_orig(:,sort(r(1:c))); 
        template = template_orig(:,sort(r(1:c)));  % sort is just for vis. later
        
%         rr = randperm(length(data(:)));
%         data(rr(1:round(i*c/10)))=1;
        rr = randperm(c);
%         data(:,r(1:round(i*c/10))) = bz_shuffleCircular(data(:,r(1:round(i*c/10)))')';
%         data(:,r(1:round(i*c/10))) = bz_shuffleCellID(data(:,r(1:round(i*c/10)))')';
        
        % normalze count to FR
        data = data./binSize;
        %% run all equal rates
        template = template_orig(:,sort(r(1:c)));  % 
        for j=1:c
            template(:,j) = minmax_norm(template(:,j))*rates_flat(r(j)); 
        end
        template(:,rr(1:round(i*c/10))) = template(:,rr(1:round(i*c/10))) +  bz_shuffleCellID(template(:,rr(1:round(i*c/10)))')';
        
        [Pr, prMax] = placeBayes(data, template', .01);
        weightedCorr(c,iter,i) = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1));
        [a b ord_template] = sort_cells(template');
        [a ord_firstSpk] = sortrows(data','descend');
        rankOrder(c,iter,i) = corr(ord_template,ord_firstSpk,'rows','complete');;
        % subplot(2,2,1)
        % plot(mean(weightedCorr,1))
        % hold on
        % plot(mean(rankOrder,1))
        % hold off

        
         %% run uniform FR distro
%         r = randperm(100);
%         data = data_orig(:,r(1:c)); template = template_orig(:,r(1:c));
%         r = randperm(length(data(:)));
%         data(r(1:round(i*c/10)))=1;
%         r = randperm(100);
        template = template_orig(:,sort(r(1:c)));  % 
        for j=1:c
            template(:,j) = minmax_norm(template(:,j))*rates_unif(r(j)); 
        end
        template(:,rr(1:round(i*c/10))) = template(:,rr(1:round(i*c/10))) +  bz_shuffleCellID(template(:,rr(1:round(i*c/10)))')';
        
        [Pr, prMax] = placeBayes(data, template', .01);
        weightedCorr_unif(c,iter,i) = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1));
        [a b ord_template] = sort_cells(template');
        [a ord_firstSpk] = sortrows(data','descend');
        rankOrder_unif(c,iter,i) = corr(ord_template,ord_firstSpk,'rows','complete');
        % subplot(2,2,2)
        % plot(mean(weightedCorr_unif,1))
        % hold on
        % plot(mean(rankOrder_unif,1))
        % hold off

        %% run log Norm FR distro
%         r = randperm(100);
%         data = data_orig(:,r(1:c)); template = template_orig(:,r(1:c));
%         r = randperm(length(data(:)));
%         data(r(1:round(i*c/10)))=1;
%         r = randperm(100);
        template = template_orig(:,sort(r(1:c)));  % 
        for j=1:c
            template(:,j) = minmax_norm(template(:,j))*rates_logN(r(j)); 
        end
        template(:,rr(1:round(i*c/10))) = template(:,rr(1:round(i*c/10))) +  bz_shuffleCellID(template(:,rr(1:round(i*c/10)))')';
        
        [Pr, prMax] = placeBayes(data, template', .01);
        weightedCorr_logN(c,iter,i) = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1));
        [a b ord_template] = sort_cells(template');
        [a ord_firstSpk] = sortrows(data','descend');
        rankOrder_logN(c,iter,i) =corr(ord_template,ord_firstSpk,'rows','complete');
        % subplot(2,2,3)
        % plot(mean(weightedCorr_logN,1))
        % hold on
        % plot(mean(rankOrder_logN,1))
        % hold off

        pause(.0001)
    end
    
%     histogram(weightedCorr(c,:));
%     histogram(weightedCorr_unif(c,:));
%     histogram(weightedCorr_logN(c,:))
    for i=1:8
    cc(c,iter,i) = corr(squeeze(rankOrder(c,iter,i:i+2)),squeeze(weightedCorr(c,iter,i:i+2)));
    c_unif(c,iter,i) = corr(squeeze(rankOrder_unif(c,iter,i:i+2)),squeeze(weightedCorr_unif(c,iter,i:i+2)));
    c_logN(c,iter,i) = corr(squeeze(rankOrder_logN(c,iter,i:i+2)),squeeze(weightedCorr_logN(c,iter,i:i+2)));
%         cc(c,iter,i) = corr(squeeze(rankOrder(c,:,i))',squeeze(weightedCorr(c,:,i))');
%         c_unif(c,iter,i) = corr(squeeze(rankOrder_unif(c,:,i))',squeeze(weightedCorr_unif(c,:,i))');
%         c_logN(c,iter,i) = corr(squeeze(rankOrder_logN(c,:,i))',squeeze(weightedCorr_logN(c,:,i))');
%     
    end
    
    subplot(4,2,1)
    imagesc(squeeze(nanmean(cc,2))); caxis([0 .75])
    subplot(4,2,2)
    imagesc(squeeze(nanmean(c_unif,2))); caxis([0 .75])
    subplot(4,2,3)
    imagesc(squeeze(nanmean(c_logN,2))); caxis([0 .75])
    subplot(4,2,4)
    plot(squeeze(nanmean(nanmean(cc))),'k')
    hold on
    plot(squeeze(nanmean(nanmean(c_unif))),'b')
    plot(squeeze(nanmean(nanmean(c_logN))),'r')
    hold off
    subplot(4,2,5)
    imagesc(squeeze(nanmean(rankOrder,2))); caxis([0 .75])
    subplot(4,2,6)
    imagesc(squeeze(nanmean(rankOrder_logN,2))); caxis([0 .75])
    subplot(4,2,7)
    imagesc(squeeze(nanmean(weightedCorr,2))); caxis([0 .75])
    subplot(4,2,8)
    imagesc(squeeze(nanmean(weightedCorr_logN,2))); caxis([0 .75])
    pause(.0001)
    end
end