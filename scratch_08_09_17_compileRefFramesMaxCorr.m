clear
d  = dir('*201*');
hpc_count = 1;
ls_count = 1;
ls_phase = []; ls_rate = [];
hpc_phase = []; hpc_rate = [];
vals = [];
corrs = [];
smoothingRange = 1:50;
%% compile data  
for ii=1:length(d)
   cd(d(ii).name)
   if exist([d(ii).name '.referenceFramesMaxCorr_fast.mat'])
   load([d(ii).name '.referenceFramesMaxCorr_fast.mat'])
   load([d(ii).name '.spikes.cellinfo.mat'])
   for cell = 1:length(mse_all_rate)
       if size(mse_all_rate{cell},1) > 3
           
       mse_rate = nanmean(mse_all_rate{cell},3);
       mse_phase = nanmean(mse_all_phase{cell},3);
       mse_phase_cos = nanmean(mse_all_phase_cos{cell},3);
       mse_phase_cos_chance = nanmean(mse_all_phase_cos_chance{cell},3);
       mse_phase_sin = nanmean(mse_all_phase_sin{cell},3);
       mse_phase_chance = nanmean(mse_all_phase_chance{cell},3);
%     mse_chance_rate = mse_all_chance_rate{cell};
%     mse_chance_phase_sin = mse_all_chance_phase_sin{cell};
%     mse_chance_phase_cos = mse_all_chance_phase_cos{cell};

%         for j=1:size(mse_rate,1)
                mse_norm_rate = (mse_rate);-nanmean(nanmean((mse_rate))); %zscore
%         end
        
%         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse_norm_rate,3))))
%         title('rate')
%                 line([6.5 6.5],[0 4000],'color','k')
%         line([10.5 10.5],[0 4000],'color','k')
%         line([18.5 18.5],[0 4000],'color','k')
%         line([19.5 19.5],[0 4000],'color','k')
        
%         subplot(2,2,4)
%         for j=1:size(mse_phase_cos,1)
                mse_norm_phase = (mse_phase_cos); %zscore
                mse_norm_phase_chance = (mse_phase_cos_chance); %zscore
%         end
%         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse_norm_phase,3))))
%         line([6.5 6.5],[0 4000],'color','k')
%         line([10.5 10.5],[0 4000],'color','k')
%         line([18.5 18.5],[0 4000],'color','k')
%         line([19.5 19.5],[0 4000],'color','k')
% %         line([6.5 6.5],[0 4000],'color','k')
%         title('phase')
%         hold off
%         
%         subplot(2,2,1)
%         % allo
%         plot(smoothingRange(1:c),((squeeze(nanmean(nanmean(mse_rate(1:c,1:6,:),2),3)))),'.k')
%         hold on
%         %route
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(1:c,7:10,:),2),3)))),'.r')
%         %cond
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(1:c,11:20,:),2),3)))),'.b')
%         %goal
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(1:c,21,:),2),3)))),'.g')
%         %ego
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_rate(1:c,22:end,:),2),3)))),'.m')
%         hold off
%         title(spikes.region{cell})
%         subplot(2,2,3)
%         % allo
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,1:6,:),2),3))))./mean(mse_chance_phase,2),'.k')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(1:c,1:6,:),2),3)))),'.k')
%         hold on
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(1:c,1:6,:),2),3)))),'.k')
%         %route
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,7:10,:),2),3))))./mean(mse_chance_phase,2),'.r')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(1:c,7:10,:),2),3)))),'.r')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(1:c,7:10,:),2),3)))),'.r')
%         %cond
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,11:20,:),2),3))))./mean(mse_chance_phase,2),'.b')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(1:c,11:20,:),2),3)))),'.b')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(1:c,11:20,:),2),3)))),'.b')
%         %goal
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,21,:),2),3))))./mean(mse_chance_phase,2),'.g')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(1:c,21,:),2),3)))),'.g')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(1:c,21,:),2),3)))),'.g')
%         %ego
% %         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase(:,22:end,:),2),3))))./mean(mse_chance_phase,2),'.m')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_cos(1:c,22:end,:),2),3)))),'.m')
%         plot(smoothingRange(1:c),((squeeze(mean(mean(mse_phase_sin(1:c,22:end,:),2),3)))),'.m')
%         hold off
        
%         imagesc(1:p,smoothingRange(1:c),(squeeze(mean(mse,3))))

%         [a b]= corr(pp(:),rr(:),'rows','complete');
%         title([cell])
%         corrs = [corrs; a sum(double(spikes.region{cell}))];
%         pause(.0001)
%        vals = [vals,a];
%        if size(pp,2) == 31
           if strcmp(spikes.region{cell},'ls')
               ls_phase(ls_count,:,:) = mse_norm_phase([1:10 end-12:end],:)';
               ls_rate(ls_count,:,:) = mse_norm_rate([1:10 end-12:end],:)';
               ls_phase_chance(ls_count,:,:) = mse_norm_phase_chance([1:10 end-12:end],:)';
               ls_count = 1+ls_count;
           elseif strcmp(spikes.region{cell},'hpc')
               hpc_phase(hpc_count,:,:) = mse_norm_phase([1:10 end-12:end],:)';
               hpc_rate(hpc_count,:,:) = mse_norm_rate([1:10 end-12:end],:)';
               hpc_count = 1+hpc_count;
           end
%        else
%            disp(['not using rec cause wrong size.. ' spikes.sessionName])
%            scratch_07_20_17_decoding_ego_allo_route_centric
%        end
       clear mse_norm_phase mse_norm_rate
       end
   end
   else
       disp(['not using ' d(ii).name])
%        scratch_07_20_17_decoding_ego_allo_route_centric
   end 
   cd('/home/david/datasets/lsDataset')
%     cd('D:\Dropbox\datasets\lsDataset')
end

% 1:6 allo 7:10 route 11 goal 12:23 ego
legend = {'x','y','z','pitch','yaw','roll','routes','goal','ego(accel)','ego(vel)'};
legendShort = {'allo','route','goal','ego'};
% figure
windRange = 21;
clf
subplot(3,1,1)
errorbar(1:23,squeeze(nanmean(nanmean(ls_phase_chance(:,windRange,:),1),2)),squeeze(nanmean(nanstd(ls_phase_chance(:,windRange,:),1),2))./sqrt(size(ls_phase,1)),'.r')
hold on
errorbar(1.5:1:23.5,squeeze(nanmean(nanmean(ls_phase(:,windRange,:),1),2)),squeeze(nanmean(nanstd(ls_phase(:,windRange,:),1),2))./sqrt(size(ls_phase,1)),'.k')
xticks([1:6 8.5 11 14.5 20.5])
set(gca,'xticklabel',legend)
title([num2str(size(ls_phase,1)) ' LS cells included'])
axis([0 24 1 1.5])
subplot(3,1,2)
errorbar(1.1,mean(squeeze(nanmean(nanmean(ls_phase_chance(:,windRange,1:3))))),(squeeze(nanstd(nanmean(ls_phase_chance(:,windRange,1:3),3))))./sqrt(size(ls_phase_chance,1)),'r')
hold on
errorbar(2.1,mean(squeeze(nanmean(nanmean(ls_phase_chance(:,windRange,7:10))))),(squeeze(nanstd(nanmean(ls_phase_chance(:,windRange,7:10),3))))./sqrt(size(ls_phase_chance,1)),'r')
errorbar(3.1,mean(squeeze(nanmean(nanmean(ls_phase_chance(:,windRange,11))))),(squeeze(nanstd(nanmean(ls_phase_chance(:,windRange,11),3))))./sqrt(size(ls_phase_chance,1)),'r')
errorbar(4.1,mean(squeeze(nanmean(nanmean(ls_phase_chance(:,windRange,12:end))))),(squeeze(nanstd(nanmean(ls_phase_chance(:,windRange,12:end),3))))./sqrt(size(ls_phase_chance,1)),'r')
% controls
errorbar(1,mean(squeeze(nanmean(nanmean(ls_phase(:,windRange,1:3))))),(squeeze(nanstd(nanmean(ls_phase(:,windRange,1:3),3))))./sqrt(size(ls_phase,1)))
hold on
errorbar(2,mean(squeeze(nanmean(nanmean(ls_phase(:,windRange,7:10))))),(squeeze(nanstd(nanmean(ls_phase(:,windRange,7:10),3))))./sqrt(size(ls_phase,1)))
errorbar(3,mean(squeeze(nanmean(nanmean(ls_phase(:,windRange,11))))),(squeeze(nanstd(nanmean(ls_phase(:,windRange,11),3))))./sqrt(size(ls_phase,1)))
errorbar(4,mean(squeeze(nanmean(nanmean(ls_phase(:,windRange,12:end))))),(squeeze(nanstd(nanmean(ls_phase(:,windRange,12:end),3))))./sqrt(size(ls_phase,1)))
xticks([1:4])
set(gca,'xticklabel',legendShort)
axis([0 5 1 1.5])

