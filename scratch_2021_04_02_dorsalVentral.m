 clear all
%% DO CCG replication here...

cd /mnt/nyuShare/Homes/dwt244/glucoseRev/dvHPC_rec/CGM37
d = dir('CGM37_*');
for i=[1:length(d)]
cd(d(i).name)
if exist([d(i).name  '.DorsalRipples.events.mat']) & exist([d(i).name  '.VentralRipples.events.mat'])
dors = load([d(i).name  '.DorsalRipples.events.mat']);
load([d(i).name  '.VentralRipples.events.mat']);
[times groups] = spikes2sorted({dors.ripples.peaks,ripples.peaks});
[ccg t] = CCG(times,groups,'binSize',.01,'duration',1,'norm','rate');
[ccg_slow t] = CCG(times,groups,'binSize',1,'duration',60*60,'norm','rate');
if size(ccg,2)>1
cc_37(i,:) = squeeze(ccg(:,1,2));
cc_37_slow(i,:) = squeeze(ccg_slow(:,1,2));
end
end
cd ..
end
imagesc(zscore(cc_37,[],2))
 plot(
cd /mnt/nyuShare/Homes/dwt244/glucoseRev/dvHPC_rec/CGM48
d = dir('CGM48_0*');
for i=[1:length(d)]
cd(d(i).name)
if exist([d(i).name  '.DorsalRipples.events.mat']) & exist([d(i).name  '.VentralRipples.events.mat'])
dors = load([d(i).name  '.DorsalRipples.events.mat']);
load([d(i).name  '.VentralRipples.events.mat']);
[times groups] = spikes2sorted({dors.ripples.peaks,ripples.peaks});
[ccg t] = CCG(times,groups,'binSize',.01,'duration',1,'norm','rate');
[ccg_slow t] = CCG(times,groups,'binSize',1,'duration',60*60,'norm','rate');
if size(ccg,2)>1
cc_48(i,:) = squeeze(ccg(:,1,2));
cc_48_slow(i,:) = squeeze(ccg_slow(:,1,2));
end
end
cd ..
end
imagesc(zscore(cc_48,[],2))
 
cd /mnt/nyuShare/Homes/dwt244/glucoseRev/dvHPC_rec/Vanessa
d = dir('Vanessa_*');
for i=[1:length(d)]
cd(d(i).name)
if exist([d(i).name  '.DorsalRipples.events.mat']) & exist([d(i).name  '.VentralRipples.events.mat'])
dors = load([d(i).name  '.DorsalRipples.events.mat']);
load([d(i).name  '.VentralRipples.events.mat']);
[times groups] = spikes2sorted({dors.ripples.peaks,ripples.peaks});
[ccg t] = CCG(times,groups,'binSize',.01,'duration',1,'norm','rate');
[ccg_slow t] = CCG(times,groups,'binSize',1,'duration',60*60,'norm','rate');
if size(ccg,2)>1
cc_van(i,:) = squeeze(ccg(:,1,2));
cc_van_slow(i,:) = squeeze(ccg_slow(:,1,2));
end
end
cd ..
end
imagesc(zscore(cc_van,[],2))




%% cross ref of dorsal/ventral 'types'?


%% glucose time
cd /media/data/Dropbox/Documents/pubs/inProgress/glucose/data/revision
r = {'CGM37.mat';'CGM48.mat';'Vanessa.mat'};
  
 for rr = 1:length(r)
    load(r{rr},'count','vcount','isig_levels','idx')
     
    clear cc ccc
    for ii=41:length(isig_levels)-40
    id = ii-40:ii+40; cc = nan(61,1);
    for i = -40:40
    if sum(~isnan(count(id)))>40
    cc(i+41) = corr(circshift(count(id),i)',[0;diff(isig_levels(id)')],'rows','complete');
    else
    cc(i+41)=nan;
    end
    end
    ccc(ii,:) = cc; clear cc;
    end
    ccc(end:length(count),:)=nan;

%     clear cc ccc
    for ii=41:length(isig_levels)-40
    id = ii-40:ii+40; cc = nan(61,1);
    for i = -40:40
    if sum(~isnan(count(id)))>40
    cc(i+41) = corr(circshift(vcount(id),i)',[0;diff(isig_levels(id)')],'rows','complete');
    else
    cc(i+41)=nan;
    end
    end
    vccc(ii,:) = cc; clear cc;
    end
    vccc(end:length(count),:)=nan;
    
    ccg_d(rr,:) = nanmedian(ccc(idx,:));
    ccg_v(rr,:) = nanmedian(vccc(idx,:));
    cc_dors{rr} = ccc(idx,:);
    cc_vent{rr} = vccc(idx,:);
    clear ccc vccc
 end

for i=1:3
    pks(i,:) = [min(ccg_v(i,40:45)) min(ccg_v(i,40:45))];
    avg(i,:) =  [mean(ccg_v(i,42:44)) mean(ccg_v(i,42:44))];
end
 
 
 %% plotting
subplot(3,2,5)
plot(nanmean(ccg_v))
hold on
plot(nanmean(ccg_v))

subplot(3,2,6)
plot(pks','k')


 cd /media/data/Dropbox/Documents/pubs/inProgress/glucose
 
 save('dorsalVentralDataset.mat','-v7.3')
