for s = 1:length(spikes.times)
    if strcmp(spikes.region{s},'ls')
        mark(s) = 1;
    else
        mark(s) = 0;
    end
end
firstHalf = [];secondHalf=[];
for lsNeuron = 62%1:sum(mark)
    for con = 1:length(unique(behavior.events.trialConditions))
    for k =1:length(spikes.times)
    if ~isempty(fields{con}{k}) &  size(firingMaps.rateMaps{con},2) > 10 
    for i=1:size(firingMaps.rateMaps{con},2)
    f = find(firingMaps.phaseMaps{con}{k}(:,2)==i);
    ff = find(firingMaps.phaseMaps{con}{k}(f,1)<fields{con}{k}{1}.COM  & fields{con}{k}{1}.start< firingMaps.phaseMaps{con}{k}(f,1));
    if ~isempty(ff) & ~isempty(firingMaps.phaseMaps{con}{lsNeuron})
    l = find(firingMaps.phaseMaps{con}{lsNeuron}(:,2)==i);
    if ~isempty(l)
    [times groups] = spikes2sorted([{firingMaps.phaseMaps{con}{lsNeuron}(l,5)},{firingMaps.phaseMaps{con}{k}(f(ff),5)}]);
    cc = CCG(times,groups,'binSize',.001,'duration',.2);
    c{con}(k,i,:) = cc(:,1,2);
    [a b] = max(squeeze(cc(:,1,2)));
    if b ~= 101
    firstHalf = [firstHalf;squeeze(cc(:,1,2))'];
    end
    end
    end
    end
    end
    end
    end
    for con = 1:length(unique(behavior.events.trialConditions))
    for k =1:length(spikes.times)
    if ~isempty(fields{con}{k}) &  size(firingMaps.rateMaps{con},2) > 10
    for i=1:size(firingMaps.rateMaps{con},2)
    f = find(firingMaps.phaseMaps{con}{k}(:,2)==i);
    ff = find(firingMaps.phaseMaps{con}{k}(f,1)>fields{con}{k}{1}.COM   & fields{con}{k}{1}.stop > firingMaps.phaseMaps{con}{k}(f,1));
    if ~isempty(ff)  & ~isempty(firingMaps.phaseMaps{con}{lsNeuron})
    l = find(firingMaps.phaseMaps{con}{lsNeuron}(:,2)==i);
    if ~isempty(l)
    [times groups] = spikes2sorted([{firingMaps.phaseMaps{con}{lsNeuron}(l,5)},{firingMaps.phaseMaps{con}{k}(f(ff),5)}]);
    cc = CCG(times,groups,'binSize',.001,'duration',.2);
    c_end{con}(k,i,:) = cc(:,1,2);
    [a b] = max(squeeze(cc(:,1,2)));
    if b ~= 101
    secondHalf = [secondHalf;squeeze(cc(:,1,2))'];
    end
    end
    end
    end
    end
    end
    end
    
%     for i = 1:length(unique(behavior.events.trialConditions))
%     for j =1:length(spikes.times)
%     if ~isempty(fields{j}{i})
%     subplot(2,2,1)
%     imagesc(squeeze(firingMaps.rateMaps{j}(i,:,:)))
%     subplot(2,2,2)
%     plot(t,smooth(squeeze(mean(c{j}(i,:,:)))))
%     hold on
%     plot(t,smooth(squeeze(mean(c_end{j}(i,:,:)))))
%     hold off;subplot(2,2,4)
%     plot(t,smooth(squeeze(mean(c{j}(i,:,:)))-(squeeze(mean(c_end{j}(i,:,:)))),3))
%     pause
%     end
%     end
end
 
% [ccg_cut t]= CCG(times,groups,'binSize',.001,'duration',.2);
% for i=1:67
% for j=1:8
% if ~isempty(fields{j}{i})
% subplot(2,2,1)
% imagesc(squeeze(firingMaps.rateMaps{j}(i,:,:)))
% subplot(2,2,2)
% plot(t,smooth(squeeze(mean(c{j}(i,:,:)))))
% hold on
% plot(t,smooth(squeeze(mean(c_end{j}(i,:,:)))))
% hold off;subplot(2,2,4)
% plot(t,smooth(squeeze(mean(c{j}(i,:,:)))-(squeeze(mean(c_end{j}(i,:,:)))),3))
% pause
% end
% end
% end