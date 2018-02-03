clear
d  = dir('*201*');
hpc_count = 1;
ls_count = 1;
ls_phase = []; ls_rate = [];
hpc_phase = []; hpc_rate = [];
vals = [];
corrs = [];
smoothingRange = 1:300:4000;
%% compile data  
for ii=length(d):-1:1
   cd(d(ii).name)
   if exist([d(ii).name '.referenceFramesMaxCorr.mat'])
   load([d(ii).name '.referenceFramesMaxCorr.mat'])
   load([d(ii).name '.spikes.cellinfo.mat'])
   for cell = 1:length(mse_all_rate)
%        alloScore = squeeze(nanmean(nanmean(mse_all_rate{cell}(:,1:2,:),3),2));
%        alloScore_chance = squeeze((nanmean(mse_all_chance_rate{cell},2)));
%        [a b] = min(alloScore./alloScore_chance);
       
       
       mse_rate = mse_all_rate{cell};
    mse_phase_cos = mse_all_phase_cos{cell};
    mse_phase_sin = mse_all_phase_sin{cell};
%     mse_chance_rate = mse_all_chance_rate{cell};
%     mse_chance_phase_sin = mse_all_chance_phase_sin{cell};
%     mse_chance_phase_cos = mse_all_chance_phase_cos{cell};
               subplot(2,2,2)
               mse_rate(isnan(mse_rate))=nanmean(mse_rate(:));
        for j=1:size(mse_rate,1)
                mse_norm_rate(j,:,:) = zscore(mse_rate(j,:,:));
        end
        p=33;c=14;
        imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse_norm_rate,3))))
        title('rate')
                line([6.5 6.5],[0 4000],'color','k')
        line([10.5 10.5],[0 4000],'color','k')
        line([18.5 18.5],[0 4000],'color','k')
        line([19.5 19.5],[0 4000],'color','k')
        
        subplot(2,2,4)
        for j=1:size(mse_phase_cos,1)
                mse_norm_phase(j,:,:) = zscore(mse_phase_cos(j,:,:));
        end
        imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse_norm_phase,3))))
        line([6.5 6.5],[0 4000],'color','k')
        line([10.5 10.5],[0 4000],'color','k')
        line([18.5 18.5],[0 4000],'color','k')
        line([19.5 19.5],[0 4000],'color','k')
%         line([6.5 6.5],[0 4000],'color','k')
        title('phase')
        hold off
        
        subplot(2,2,1)
        % allo
        plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_rate(:,1:6,:),2),3)))),'.k')
        hold on
        %route
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,7:10,:),2),3)))),'.r')
        %cond
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,11:20,:),2),3)))),'.b')
        %goal
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,21,:),2),3)))),'.g')
        %ego
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(:,22:end,:),2),3)))),'.m')
        hold off
        title(spikes.region{cell})
        subplot(2,2,3)
        % allo
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,1:6,:),2),3))))./mean(mse_chance_phase,2),'.k')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(:,1:6,:),2),3)))),'.k')
        hold on
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(:,1:6,:),2),3)))),'.k')
        %route
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,7:10,:),2),3))))./mean(mse_chance_phase,2),'.r')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(:,7:10,:),2),3)))),'.r')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(:,7:10,:),2),3)))),'.r')
        %cond
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,11:20,:),2),3))))./mean(mse_chance_phase,2),'.b')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(:,11:20,:),2),3)))),'.b')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(:,11:20,:),2),3)))),'.b')
        %goal
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,21,:),2),3))))./mean(mse_chance_phase,2),'.g')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(:,21,:),2),3)))),'.g')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(:,21,:),2),3)))),'.g')
        %ego
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,22:end,:),2),3))))./mean(mse_chance_phase,2),'.m')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(:,22:end,:),2),3)))),'.m')
        plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(:,22:end,:),2),3)))),'.m')
        hold off
        
%         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse,3))))
pp = (squeeze(nanmean(mse_norm_phase,3)));
rr = (squeeze(nanmean(mse_norm_rate,3)));
        [a b]= corr(pp(:),rr(:),'rows','complete');
        title([cell a])
        corrs = [corrs; a sum(double(spikes.region{cell}))];
        pause(.0001)
       vals = [vals,a];
%        if size(pp,2) == 31
           if strcmp(spikes.region{cell},'ls')
               ls_phase(ls_count,:,:) = pp(:,[1:10 end-13:end]);
               ls_rate(ls_count,:,:) = rr(:,[1:10 end-13:end]);
               
               ls_count = 1+ls_count;
           elseif strcmp(spikes.region{cell},'hpc')
               hpc_phase(hpc_count,:,:) = pp(:,[1:10 end-13:end]);
               hpc_rate(hpc_count,:,:) = rr(:,[1:10 end-13:end]);
               
               hpc_count = 1+hpc_count;
           end
%        else
%            disp(['not using rec cause wrong size.. ' spikes.sessionName])
%            scratch_07_20_17_decoding_ego_allo_route_centric
%        end
       clear mse_norm_phase mse_norm_rate
   end
   else
%        scratch_07_20_17_decoding_ego_allo_route_centric
   end 
   cd('/home/david/datasets/lsDataset')
end
