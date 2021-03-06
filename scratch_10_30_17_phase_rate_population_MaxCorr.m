clf
% cd D:\Dropbox\datasets\lsDataset
cd /home/david/datasets/lsDataset
clear all
warning off
ii=1
for tau = 1:10
%     figure(tau)
    
    d = dir('*201*');
    for ii=length(d):-1:1
        cd(d(ii).name)
        if exist([d(ii).name '.firingMaps.cellinfo.mat'])% & exist([d(ii).name '.placeFields.10_pctThresh.mat'])
        sessionInfo = bz_getSessionInfo;
        load([sessionInfo.FileName '.behavior.mat'])
        load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
%         load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])

        spikes = bz_GetSpikes;
        mark=[];
        for i=1:length(spikes.times)
            if  strcmp(spikes.region{i},'ls') %strcmp(spikes.region{i},'hpc') || strcmp(spikes.region{i},'ca3') ||strcmp(spikes.region{i},'ca1') % strcmp(spikes.region{i},'ls') %
                mark(i) = 1;
            else
                mark(i) = 0;
            end
        end
        [phaseMaps] = bz_phaseMap2Bins(firingMaps.phaseMaps,firingMaps.rateMaps,behavior);
        [firingMaps] = bz_firingMap1D(spikes,behavior,tau);
        % phase smoothing here
        for i=1:length(behavior.events.conditionType)
        for cell=1:length(spikes.times)
        for k=1:size(phaseMaps{i},2)
        pp{i}(cell,k,:) = circ_smoothTS(squeeze(phaseMaps{i}(cell,k,:))',tau,'method','mean');
        end
        end
        end

        for i=1:length(behavior.events.conditionType)
        for cell=1:length(spikes.times)
            p{i}(cell,:) = reshape(squeeze(pp{i}(cell,:,:))',1,size(pp{i},3)*size(pp{i},2));

            % rate smoothing here
            for trial = 1:size(firingMaps.rateMaps{i},2)
            smoothed_rates(trial,:) = smooth(squeeze(firingMaps.rateMaps_box{i}(cell,trial,:))',tau);
            end
            r{i}(cell,:) = reshape(smoothed_rates',1,size(pp{i},3)*size(firingMaps.rateMaps{i},2));
            clear smoothed_rates
        end
        t{i} = (repmat([1:size(pp{i},3)]',size(pp{i},2),1));
        end

        for cond = 1:length(behavior.events.conditionType)
            for cell =1:length(spikes.times)
            mark_pf(cell) = 1;%~isempty(fields{cond}{cell}); 
            end
            if length(intersect(find(mark==1),find(mark_pf==1))) >= 5 && size(firingMaps.rateMaps{cond},2) > 12

            pos = t{cond};
            rate = r{cond}(intersect(find(mark==1),find(mark_pf==1)),:);
            phase = p{cond}(intersect(find(mark==1),find(mark_pf==1)),:);  

            for iter = 1:5
                rr = randperm(length(pos));
                r_train = rr(1:prctile(1:length(pos),60));
                r_test = rr(prctile(1:length(pos),60):length(pos));

                cl = max_correlation_coefficient_CL;
                cl = train(cl,rate(:,r_train),pos(r_train));
                yfit_rate = test(cl,rate(:,r_test));
                mse_rate{ii}(cond,iter) = mean((yfit_rate-pos(r_test)).^2);

                cl = max_correlation_coefficient_CL;
                cl = train(cl,phase(:,r_train),pos(r_train));
                yfit_phase = test(cl,phase(:,r_test));
                mse_phase{ii}(cond,iter) = mean((yfit_phase-pos(r_test)).^2);

                cl = max_correlation_coefficient_CL;
                phase_c = cos(phase);
                cl = train(cl,phase_c(:,r_train),pos(r_train));
                yfit_phase_c = test(cl,phase_c(:,r_test));
                mse_phase_c{ii}(cond,iter) = mean((yfit_phase_c-pos(r_test)).^2);

                cl = max_correlation_coefficient_CL;
                phase_s = cos(phase);
                cl = train(cl,phase_s(:,r_train),pos(r_train));
                yfit_phase_s = test(cl,phase_s(:,r_test));
                mse_phase_s{ii}(cond,iter) = mean((yfit_phase_s-pos(r_test)).^2);

                cl = max_correlation_coefficient_CL;
                cl = train(cl,[phase_s(:,r_train); phase_c(:,r_train)],pos(r_train));
                yfit_phase_b = test(cl,[phase_s(:,r_test); phase_c(:,r_test)]);
                mse_phase_b{ii}(cond,iter) = mean((yfit_phase_b-pos(r_test)).^2);

            end
%                 subplot(3,2,1)
%                 scatter(mean(mse_phase{ii}(cond,:),2),mean(mse_rate{ii}(cond,:),2),'.k')
%                 hold on
%                 set(gca,'xscale','log')
%                 set(gca,'yscale','log')
% 
%                 subplot(3,2,2)
%                 scatter(mean(mse_phase_c{ii}(cond,:),2),mean(mse_rate{ii}(cond,:),2),'.k')
%                 title('cos')
%                 hold on
%                 set(gca,'xscale','log')
%                 set(gca,'yscale','log')
% 
%                 subplot(3,2,3)
%                 scatter(mean(mse_phase_s{ii}(cond,:),2),mean(mse_rate{ii}(cond,:),2),'.k')
%                 title('sin')
%                 hold on
%                 set(gca,'xscale','log')
%                 set(gca,'yscale','log')
% 
%                 subplot(3,2,4)
%                 scatter(mean(mse_phase_b{ii}(cond,:),2),mean(mse_rate{ii}(cond,:),2),'.k')
%                 title('both')
%                 hold on
%                 set(gca,'xscale','log')
%                 set(gca,'yscale','log')
% 
%                 subplot(3,2,5)
%                 numCells{ii}(cond) = length(intersect(find(mark==1),find(mark_pf==1)));
%                 hold on
%                 scatter(numCells{ii}(cond),log(mean(mse_rate{ii}(cond,:),2))-log(mean(mse_phase_c{ii}(cond,:),2)),'.k')
% 
%                 subplot(3,2,6)
%                 scatter(numCells{ii}(cond),mean(mse_rate{ii}(cond,:),2),'.k')
%                 scatter(numCells{ii}(cond),mean(mse_phase{ii}(cond,:),2),'.g')
%                 hold on
%                 set(gca,'yscale','log')
                scatter(tau,log(mean(mse_rate{ii}(cond,:),2))-log(mean(mse_phase_c{ii}(cond,:),2)),'.k')
                hold on
                pause(.01)
            end
        end
        end
        clear r p t pp pos rate phase mark
%         cd D:\Dropbox\datasets\lsDataset
        cd /home/david/datasets/lsDataset
    end
end