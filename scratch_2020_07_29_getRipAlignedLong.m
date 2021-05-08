cd('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data\stim')
d = dir('cgm*')

for thresh = 40
for i=1:8
load(d(i).name,'cc','absTime','isa','ripTime','idx','ccc','isig_levels','count','stim*','zt','responseMag');
stimAlignedL = nan(length(stimTime),thresh*2 +1);
parfor ii=1:length(stimTime)
[a b]=min(abs(stimTime(ii)-absTime));
aa(ii)=a;
if b>thresh & b<length(absTime)-thresh
stimAlignedL(ii,:)=isig_levels(b-thresh:b+thresh);
end
end
id = find((absTime(idx(1))<stimTime & absTime(idx(end))>stimTime));
c(i,:)=nanmean(diff(stimAlignedCGM(id,:)')');
cl(i,:)=nanmean(diff(stimAlignedL(id,:)')');
stimsAligned{i} = diff(stimAlignedL(id,:)')';
respsAligned{i} = responseMag(id,:);
clear stimAligned* responseMag;
i
end
end
whos resps* stims*
% 
% for thresh = 1:10000
% for k=1:2
% for j =1:8
% for i=1:80
% [a(j,i,k) b(j,i,k)] = corr((respsAligned{j}(:,k)-movmean(respsAligned{j}(:,k),thresh)),nanmean(stimsAligned{j}(:,i),2),'rows','complete');
% end
% end
% end
% for j =1:8
% az(thresh,j,:)=nanZscore(a(j,:,ss(j)));
% aa(thresh,j,:)=(a(j,:,ss(j)));
% bb(thresh,j,:)=b(j,:,ss(j));
% end; aAll(thresh,j,:,:)= a(j,:,:);
% imagesc(squeeze(mean(aa,2)))
% pause(.1)
% end

cd('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data')

d = dir('_*mat');
for thresh = 40
for i=1:8
load(d(i).name,'cc','absTime','isa','rip*','idx','ccc','isig_levels','count','stim*','zt','responseMag');
ripAlignedL = nan(length(ripTime),thresh*2 +1);
parfor ii=1:length(ripTime)
[a b]=min(abs(ripTime(ii)-absTime));
aa(ii)=a;
if b>thresh & b<length(absTime)-thresh
ripAlignedL(ii,:)=isig_levels(b-thresh:b+thresh);
end
end
id = find((absTime(idx(1))<ripTime & absTime(idx(end))>ripTime));
c_nostim(i,:)=nanmean(diff(rippleAlignedCGM(id,:)')');
cl_nostim(i,:)=nanmean(diff(ripAlignedL(id,:)')');
ripAligned{i} = (ripAlignedL(id,:)')';
% respsAligned{i} = responseMag(id,:);
clear ripAlignedL responseMag;
i
end
end
whos resps* stims*