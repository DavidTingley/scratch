cd /home/david/datasets/lsDataset
clear
d  = dir('*201*');
count = 1;
 
%% compile data

for ii=1:length(d)
   cd(d(ii).name)
   if exist([d(ii).name '.olypherInfo.cellinfo.mat'])
   load([d(ii).name '.olypherInfo.cellinfo.mat'],'olypherInfo')  
   load([d(ii).name '.firingMaps.cellinfo.mat'],'firingMaps') 
   load([d(ii).name '.behavior.mat']) 
   nBins = round(max(behavior.events.trials{1}.mapping));
   if exist([d(ii).name '.placeFields.mat'])
   load([d(ii).name '.placeFields.mat'])
   for cell=1:length(olypherInfo.results)
       for cond = 1:length(unique(olypherInfo.results{cell}.condition))
           nTrials = size(firingMaps.rateMaps{cond},2);
%            if ~isempty(fields{cond}{cell})
            if strcmp(olypherInfo.region{cell},'hpc')
                
            disc = unique(olypherInfo.results{cell}.discBins);
            for dd = 5%:length(disc)
                
            rows = find(olypherInfo.results{cell}.condition==cond);
            cols = find(olypherInfo.results{cell}.discBins==disc(dd));
            rows = intersect(rows,cols);
            
            if length(rows)>40
                % subplot(2,2,1);
                % imagesc(squeeze(rateMap{cond}(cell,:,:)));
%                 subplot(2,2,2);
%                 imagesc(squeeze(phaseMap{cond}(cell,:,:)));
                % subplot(2,2,4)
                % if ~isempty(phaseMap{cond}{cell})
                % scatter(phaseMap{cond}{cell}(:,1),phaseMap{cond}{cell}(:,end),'.k'); hold on; 
                % scatter(phaseMap{cond}{cell}(:,1),phaseMap{cond}{cell}(:,end)+2*pi,'.k'); hold off
                % end
%                 subplot(2,2,3);
% 
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.ratePeakInfo(rows),'r')
%                 hold on
%                 plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.phasePeakInfo(rows),'g')
%                 hold off

                phaseInfoScore(count,dd,:) = olypherInfo.results{cell}.phasePeakInfo(rows(1:40))./nTrials;
                rateInfoScore(count,dd,:) = olypherInfo.results{cell}.ratePeakInfo(rows(1:40))./nTrials;
                 if ~isempty(fields{cond}{cell})
                     hasField(count) = 1;
                     bins(count) = nBins;
                 else
                     hasField(count) = 0;
                     bins(count) = nBins;
                 end
          
                count = 1+count
            pause(.01)
            end
            end
           end
%            end
       end
   end
   end
   end
   cd /home/david/datasets/lsDataset
end

% subplot(2,2,1)
% boundedline(1:size(phaseInfoScore,2),smooth(nanmean(phaseInfoScore),3),nanstd(phaseInfoScore)./nanmean(phaseInfoScore),'g')
% boundedline(1:size(rateInfoScore,2),smooth(nanmean(rateInfoScore),3),nanstd(rateInfoScore)./nanmean(rateInfoScore),'r')
% 
% subplot(2,2,2)
% plot(mean(phaseInfoScore)./mean(rateInfoScore),'k')

