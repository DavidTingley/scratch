cd /home/david/datasets/lsDataset
% cd('D:\Dropbox\datasets\lsDataset')
% % clear
d  = dir('*201*');
count = 1;
 
%% compile data

for ii=1:length(d)
   cd(d(ii).name)
   if exist([d(ii).name '.olypherInfo_w_disc.cellinfo.mat'])
   load([d(ii).name '.olypherInfo_w_disc.cellinfo.mat'],'olypherInfo')  
   load([d(ii).name '.firingMaps.cellinfo.mat'],'firingMaps') 
   load([d(ii).name '.behavior.mat']) 
   sessionInfo = bz_getSessionInfo;
   nBins = length((behavior.events.trials{1}.mapping));
   if exist([d(ii).name '.placeFields.20_pctThresh.mat'])
   load([d(ii).name '.placeFields.20_pctThresh.mat'])
   conditions = length(unique(behavior.events.trialConditions));
   
   for cell=1:size(firingMaps.rateMaps{1},1)
       for cond = 1:conditions
           if sum(behavior.events.trialConditions==cond) >= 10
           nTrials = size(firingMaps.rateMaps{cond},2);
%            if ~isempty(fields{cond}{cell})
            if strcmp(olypherInfo.region{cell},'ls')% | strcmp(olypherInfo.region{cell},'ca1') | strcmp(olypherInfo.region{cell},'ca3')
                %  strcmp(olypherInfo.region{cell},'ls')%
            disc = unique(olypherInfo.results{cell}.discBins);
            if length(disc) >= 1
            for dd = 1%1:length(disc)
                
            rows = find(olypherInfo.results{cell}.condition==cond);
            cols = find(olypherInfo.results{cell}.discBins==disc(dd));
            rows = intersect(rows,cols);
            
            if length(rows)>=50
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
                if ~isempty(sessionInfo.ca3)
                    region(count) = 3;
                else
                    region(count) = 1;
                end
                phaseInfoScore(count,dd,:) = olypherInfo.results{cell}.phasePeakInfo(rows(1:50))./nTrials;
                rateInfoScore(count,dd,:) = olypherInfo.results{cell}.ratePeakInfo(rows(1:50))./nTrials;
                 if ~isempty(fields{cond}{cell})
                     hasField(count) = 1;
                     bins(count) = nBins;
                 else
                     hasField(count) = 0;
                     bins(count) = nBins;
                 end
           recording(count) =  ii;
                
            pause(.01)
%             if dd == length(disc)
            count = 1+count;
%             end
            else%if cell == 1 & cond == 1
                warning(['not using ' d(ii).name])
            end
            end
            else
                error('missing something..')
            end
           end
           end
       end
   end
   end
   else
       warning(['missed ' d(ii).name])
   end
   cd /home/david/datasets/lsDataset
   clear olypherInfo
   ii
% cd('D:\Dropbox\datasets\lsDataset')
end

% subplot(2,2,1)
% boundedline(1:size(phaseInfoScore,2),smooth(nanmean(phaseInfoScore),3),nanstd(phaseInfoScore)./nanmean(phaseInfoScore),'g')
% boundedline(1:size(rateInfoScore,2),smooth(nanmean(rateInfoScore),3),nanstd(rateInfoScore)./nanmean(rateInfoScore),'r')
% 
% subplot(2,2,2)
% plot(mean(phaseInfoScore)./mean(rateInfoScore),'k')

