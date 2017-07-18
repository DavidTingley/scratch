% compiling cell assembly data across regions and recordings...
d = dir('*201*');
clf
pairCount = 0;
h=[];
z=[];
ii=[];
rate=[];
for i=1:length(d)
cd(d(i).name)
if exist('assembliesCrossRegionData_w_theta_sin_cos_coord_vel.mat')
load('assembliesCrossRegionData_w_theta_sin_cos_coord_vel.mat')
% if exist(['assembliesCrossRegionData.mat'])
% load(['assembliesCrossRegionData.mat'])
p=[];
if ~isempty(pairs)
% figure(1);
for c = 1:length(dev)
   for pair = 1:size(dev{c},2)
    [a b] =  min(dev{c}(:,pair));
    [aa bb] = min(mean(devControl{c}(:,pair,:),3));
    imp = (a-mean(dev{c}(:,pair))) ./ (aa - mean(mean(devControl{c}(:,pair,:),3)));
    imp2 = a ./ max(mean(mean(devControl{c}(:,pair,:),3)));
    zerolag = (min(dev{c}(1:6,pair)) - mean(dev{c}(:,pair))) ./ (aa - mean(mean(devControl{c}(1,pair,:),3)));
    if zerolag < 1 
        zerolag = 1;
    end
    
    if imp > 7 & b > 7 & zerolag < 1.2 & mean(dev{c}(:,pair))>150
       subplot(2,2,1)
       scatter(b,imp,'.k')
       hold on
       subplot(2,2,2)
       histogram(ii,1:200)
       p = [p; pairs(pair,:)];
       h = [h; imp];
       ii = [ii;b];
       z=[z;zerolag];
       rate=[rate;mean(dev{c}(:,pair))];
       line([median(ii) median(ii)],[0 25],'color','r')
       line([mean(ii) mean(ii)],[0 25],'color','g')
       ylabel('improvement (min-mean) ./ (minc-meanc)')
       xlabel('time, ms')
       subplot(2,2,3)
       scatter(b,imp2,'.k')
       hold on
       subplot(2,2,4)
       plot(dev{c}(:,pair));
       hold on
       plot(mean(devControl{c}(:,pair,:),3));
       title([imp zerolag ])
       hold off
       pause(.001)
    end
    pairCount = 1 + pairCount;
   end    
end
end
end

% let's go hunting for assemblies with phase precession clouds...
% spikes = bz_GetSpikes;
% load([d(i).name '.behavior.mat'])
% load([d(i).name '.sessionInfo.mat'])
% 
% for ind=1:length(spikes.times)
%     if strcmp(spikes.region{ind},'hpc')
%         hpc(ind)=1;
%     else
%         hpc(ind)=0;
%     end
% end
% 
% theta_ind = intersect(find(hpc),sessionInfo.thetaChans);
% lfp = bz_GetLFP(theta_ind);
% 
% [rateMap countMap occuMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,4);
% 
% for ind =1:size(p)
%     figure(2)
%     for j=1:8
%     subplot(4,2,j)
%         if ~isempty(phaseMap{j})
%         if ~isempty(phaseMap{j}{p(ind,1)})
%             subplot(4,2,j)
%             scatter(phaseMap{j}{p(ind,1)}(:,1),phaseMap{j}{p(ind,1)}(:,end)+2*pi,'.k');
%             hold on
%             scatter(phaseMap{j}{p(ind,1)}(:,1),phaseMap{j}{p(ind,1)}(:,end),'.k');
%             hold off
%         end
%         end
%     end
%     figure(3)
%     for j=1:8
%         if ~isempty(phaseMap{j})
%         if ~isempty(phaseMap{j}{p(ind,2)})
%             subplot(4,2,j)
%             scatter(phaseMap{j}{p(ind,2)}(:,1),phaseMap{j}{p(ind,2)}(:,end)+2*pi,'.k');
%             hold on
%             scatter(phaseMap{j}{p(ind,2)}(:,1),phaseMap{j}{p(ind,2)}(:,end),'.k');
%             hold off
%         end
%         end
%     end
%     pause
% end




cd /home/david/datasets/lsDataset

end
