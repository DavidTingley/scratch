d = dir('*201*');
count = 1;


for i=1:length(d)
    cd(d(i).name)
    load([d(i).name '.behavior.mat'])
    if strcmp(behavior.description,'central alternation')
        totalDist = 270;
    elseif strcmp(behavior.description,'linear')
        totalDist = 200;
    elseif strcmp(behavior.description,'wheel alternation')
        totalDist = 320;
    elseif strcmp(behavior.description,'both alternation')
        totalDist = 300;
    end
    load([d(i).name '.firingMaps.cellinfo.mat'])
    if exist([d(i).name '.olypherInfo_w_disc.cellinfo.mat'])
    load([d(i).name '.olypherInfo_w_disc.cellinfo.mat'])
    load([d(i).name '.placeFields.mat'])
    nBins = size(firingMaps.rateMaps{1},3);
    for cond = 1:length(firingMaps.rateMaps)
        for cell = 1:length(olypherInfo.UID)
           if strcmp(olypherInfo.region{cell},'hpc')
           if ~isempty(fields{cond}{cell})
               rows = find(olypherInfo.results{cell}.condition==cond);
               cols = find(olypherInfo.results{cell}.discBins==10);
               rows = intersect(rows,cols);
               if ~isempty(rows)
                   
               [rate_max(count) rate_loc(count)] = max(smooth(olypherInfo.results{cell}.ratePeakInfo(rows),3));
               [phase_max(count) phase_loc(count)] = max(smooth(olypherInfo.results{cell}.phasePeakInfo(rows),3));
               
               width(count) = fields{cond}{cell}{1}.width;
               bins(count) = nBins;
               dists(count) = totalDist;
               nTrials(count) = size(firingMaps.rateMaps{cond},2);

               pause(.01)
               count = 1+count;
               end
           end           
           end       
        end
    end
    end
    
    
    cd /home/david/datasets/lsDataset/
end

subplot(3,2,1)
scatter(width./bins,rate_loc./bins,'.r')
xlabel('pf width')
ylabel('optimal smooth wind')

title('loc')
subplot(3,2,2)
scatter(width./bins,phase_loc./bins,'.g')

title('loc')
subplot(3,2,3)
scatter(width./bins,rate_max./nTrials,'.r')

title('max')
subplot(3,2,4)
scatter(width./bins,phase_max./nTrials,'.g')

title('max')
subplot(3,2,5)
scatter(phase_loc./bins.*dists,phase_max./nTrials,'.g')

title('max')
subplot(3,2,6)
scatter(rate_loc./bins.*dists,rate_max./nTrials,'.r')