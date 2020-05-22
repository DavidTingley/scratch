% need to re-do ratemapping here...
Tau_smooth = 30;

for i=1:length(rateMap)
    % set up iteration loop here..
    for iter = 1:5
        trial_ord = randperm(size(rateMap{i},2));
        rates = [];
        counts = [];
        phases_train = [];
        for j=1:size(rateMap{i},2)-1
            rates=[rates,squeeze(rateMap{i}(:,trial_ord(j),:))];
            c = squeeze(countMap{i}(:,trial_ord(j),:));
            for t = 1:size(c,1)
               c(t,:) = round(smooth(c(t,:),Tau_smooth)*Tau_smooth); 
            end
            counts = [counts,c];
            
            phases_trial = zeros(length(phaseMap{i}),(size(rateMap{1},3)));
            for cell = 1:length(phaseMap{i})
            for k = 1: size(rateMap{1},3)
                if ~isempty(phaseMap{i}{cell})
                f = find(phaseMap{i}{cell}(:,1)==k); 
                ff = find(phaseMap{i}{cell}(:,2)==j);
                if ~isempty(intersect(f,ff))
                phases_trial(cell,k) = circ_median(phaseMap{i}{cell}(intersect(f,ff),end));
                end
                end
            end
%             phases_trial(cell,:) = circ_smoothTS(phases_trial(cell,:),Tau_smooth,'exclude',0,'method','mean');
            end
            phases_train = [phases_train, phases_trial];
        end
        rates_train = rates;
        counts_train = counts;
        rates_test = squeeze(rateMap{i}(:,trial_ord(j+1),:));
        counts_test = squeeze(countMap{i}(:,trial_ord(j+1),:));
        phases = zeros(size(rateMap{1},3),1);
        phases_test = zeros(length(phaseMap{i}),(size(rateMap{1},3)));
        for cell = 1:length(phaseMap{i})
        for k = 1: size(rateMap{1},3)
            if ~isempty(phaseMap{i}{cell})
            f = find(phaseMap{i}{cell}(:,1)==k); 
            ff = find(phaseMap{i}{cell}(:,2)==j+1);
            if ~isempty(intersect(f,ff))
            phases_test(cell,k) = circ_mean(phaseMap{i}{cell}(intersect(f,ff),end));
            end
            end
        end
        end
        pos_train = repmat(1:size(rateMap{1},3),1,j);
        pos_test = 1:size(rateMap{1},3);
        
        phases_test(phases_test==0)=nan;
        phases_train(phases_train==0)=nan;
        [phases_test, EDGES] = discretize(phases_test,-pi:.5:pi);
        [phases_train, EDGES] = discretize(phases_train,-pi:.5:pi);
        phases_test(isnan(phases_test))=0;
        phases_train(isnan(phases_train))=0;
        
%         %% GLM here
%         [b dev stats] = glmfit(counts_train',pos_train','normal');
%         yfit = glmval(b,counts_test','identity');
%         mse_count_glm(i,iter) = mean(abs(yfit-pos_test')); 
% 
%         [b dev stats] = glmfit(phases_train',pos_train','normal');
%         yfit = glmval(b,phases_test','identity');
%         mse_phase_glm(i,iter) = mean(abs(yfit-pos_test')); 
%         
%         %% Bayesion here
%         cl = poisson_naive_bayes_CL;
%         cl = train(cl,(counts_train),pos_train);
%         [predicted_labels decision_values] = cl.test((phases_test));
%         mse_count_bayes(i,iter) = mean(abs(predicted_labels-pos_test)); 
%         
%         cl = poisson_naive_bayes_CL;
%         cl = train(cl,(phases_train),pos_train);
%         [predicted_labels decision_values] = cl.test((phases_test));
%         mse_phase_bayes(i,iter) = mean(abs(predicted_labels-pos_test)); 
%         
%         %% SVM
% %         cl = libsvm_CL;
% %         cl = train(cl,round(rates_train*.1),pos_train);
% %         [predicted_labels decision_values] = cl.test(round(rates_test*.1));
% %         mse_svm(i,iter) = mean(abs(predicted_labels-pos_test)); 
%         
%         %% MAX CORRELATION
%         cl = max_correlation_coefficient_CL;
%         cl = train(cl,counts_train,pos_train);
%         [predicted_labels decision_values] = cl.test(counts_test);
%         mse_count_maxCorr(i,iter) = mean(abs(predicted_labels-pos_test)); 
%         
%         cl = max_correlation_coefficient_CL;
%         cl = train(cl,phases_train,pos_train);
%         [predicted_labels decision_values] = cl.test(phases_test);
%         mse_phase_maxCorr(i,iter) = mean(abs(predicted_labels-pos_test));
% 
%         %% NNET here
%         
%         %% nSTAT methods here
%         
%     end
%     
%    subplot(2,2,3)
%    plot(median(mse_count_glm'),'.g')
%    hold on
%    plot(median(mse_count_bayes'),'.r')
%    plot(median(mse_count_maxCorr'),'.k')
%    title('counts')
% %    plot(median(mse_count_'),'.g')
% %    plot(median(mse_count_glm'),'.g')
%    subplot(2,2,4) 
%    plot(median(mse_phase_glm'),'.g')
%    hold on
%    plot(median(mse_phase_bayes'),'.r')
%    plot(median(mse_phase_maxCorr'),'.k')
%    title('phases')
%    pause(.01)
   



f = find(spikes.shankID<5);
phases_train_ls = phases_train(f,:);
rates_train_ls = rates_train(f,:);

f = find(spikes.shankID>4);
phases_train_hpc = phases_train(f,:);
rates_train_hpc = rates_train(f,:);

for j=1:size(phases_train_ls,1)
    for len = 1:200
        for ii =1:size(phases_train_hpc,1)
        phases_train_hpc_smooth(ii,:) = circ_smoothTS(phases_train_hpc(ii,:),len,'method','mean');
        end
        [b dev stats] = glmfit(phases_train_hpc_smooth',rates_train_ls(j,:)','normal');
        yfit = glmval(b,phases_train_hpc_smooth','identity');
        mse_phase(j,len) = mean((yfit-rates_train_ls(j,:)').^2);
    end
    plot(mse_phase(j,:),'r')
    hold on
    
    for len = 1:200
        for ii =1:size(rates_train_hpc,1)
        rates_train_hpc_smooth(ii,:) = smooth(rates_train_hpc(ii,:),len);
        end
        [b dev stats] = glmfit(rates_train_hpc_smooth',rates_train_ls(j,:)','normal');
        yfit = glmval(b,rates_train_hpc_smooth','identity');
        mse_rate(j,len) = mean((yfit-rates_train_ls(j,:)').^2);
    end
    plot(mse_rate(j,:),'k')
    hold off
    pause(.01)
end


    end
end

