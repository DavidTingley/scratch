cgm_files = dir('_*mat');
names = {};
for ff=1:length(cgm_files)
nostim{ff} = load(cgm_files(ff).name,'ripTime','rippleAlignedCGM','ccc','idx','isa','absTime','zt','idx','count','isig_levels');
cc_nostim(ff,:) = ccgBinned(nostim{ff}.count(nostim{ff}.idx),[0,diff(nostim{ff}.isig_levels(nostim{ff}.idx))],40);
nostim{ff}.ccc(end:length(nostim{ff}.count),:)=nan;
cc_avg_nostim(ff,:) = nanmean(nostim{ff}.ccc(nostim{ff}.idx,:));
id = find(nostim{ff}.absTime(nostim{ff}.idx(1)) < nostim{ff}.ripTime & ...
nostim{ff}.absTime(nostim{ff}.idx(end)) > nostim{ff}.ripTime);
isig = [0,diff(nostim{ff}.isig_levels)];
idx = intersect(find(~isnan(isig)),nostim{ff}.idx);
ae(ff) = approximateEntropy(isig(idx));
acg(ff,:) = ccgBinned(isig(idx),isig(idx),40);

end
cd revision
d = dir('*mat');
for j=1:length(d)
load(d(j).name,'count','isig*','idx','rippleAlignedCGM','ripTime','isa','absTime')
if strcmp(d(j).name,'CGM47.mat') | strcmp(d(j).name,'CGM50.mat')
isig = nanZscore(isa);
elseif strcmp(d(j).name,'ros.mat')
isig = isa;
else
isig = [0,diff(isig_levels)];
end

isigs{j} = isig;
if isvarname('idx') & ~isempty(idx)
idx = intersect(find(~isnan(isig)),idx);
ae(ff+j) = approximateEntropy(isig(idx));
acg(ff+j,:) = ccgBinned(isig(idx),isig(idx),40);
end
clear cc ccc
for ii=41:length(isig_levels)-40
    id = ii-40:ii+40; cc = nan(61,1);
    for i = -40:40
    if sum(~isnan(count(id)))>20
    cc(i+41) = corr(circshift(count(id),i)',[isig(id)'],'rows','complete');
    else
    cc(i+41)=nan;
    end
    end
    ccc(ii,:) = cc; clear cc;
end
whos ccc
ccc(end:length(count),:)=nan;

if isvarname('idx')
corr_th(j,:) = nanmedian(ccc(idx,:));
end

clear idx
end
cc = [corr_th; cc_avg_nostim];
whos nostim isigs cc