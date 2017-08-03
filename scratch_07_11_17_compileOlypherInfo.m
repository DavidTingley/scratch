cd /home/david/datasets/lsDataset
clear
d  = dir('*201*');
count = 1;
 
%% compile data

for ii=1:length(d)
   cd(d(ii).name)
   load([d(ii).name '.olypherInfo.cellinfo.mat'],'olypherInfo')   
   
   for cell=1:length(olypherInfo.results)
       for cond = 1:length(unique(olypherInfo.results{cell}.condition))
                subplot(2,2,1);
                imagesc(squeeze(rateMap{cond}(cell,:,:)));
%                 subplot(2,2,2);
%                 imagesc(squeeze(phaseMap{cond}(cell,:,:)));
                subplot(2,2,4)
                if ~isempty(phaseMap{cond}{cell})
                scatter(phaseMap{cond}{cell}(:,1),phaseMap{cond}{cell}(:,end),'.k'); hold on; 
                scatter(phaseMap{cond}{cell}(:,1),phaseMap{cond}{cell}(:,end)+2*pi,'.k'); hold off
                end
                subplot(2,2,3);
                rows = find(olypherInfo.results{cell}.condition==cond);
                plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.rateTotalInfo(rows),'r')
                hold on
                plot(olypherInfo.results{cell}.smoothing(rows),olypherInfo.results{cell}.phaseTotalInfo(rows),'g')
                hold off
                pause
       end
   end
   
   cd /home/david/datasets/lsDataset
   
end