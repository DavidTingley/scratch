
for i=1:782
for j=1
for k=1:150
sprc(i,j,k) = mean(success_power_rate_corrs{i,j,k});
sprc_shuffle(i,j,k,:) = mean(success_power_rate_corrs_shuffle{i,1,k}');
[a pvals(i,j,k)] = kstest2(success_power_rate_corrs{i,j,k},(success_power_rate_corrs_shuffle{i,j,k}));
r = ceil(rand*10);
[a pvals_control(i,j,k)] = kstest2(success_power_rate_corrs_shuffle{i,j,k}(r),(success_power_rate_corrs_shuffle{i,j,k}([1:r-1 r+1:end])));
end
end
end

clear r rr
for i=1:782
r(i,:) = mean(trials_all_succ{i});
end

for i=1:782
[rr(i,:)] = minmax_norm(r(i,:));
end
[a b o] = sort_cells(rr,r,1);

s = squeeze(sprc(o,1,:));
sc = squeeze((nanmean(sprc_shuffle(o,:,:,:),4)));
sc_std = squeeze((nanstd(sprc_shuffle(o,:,:,:),[],4)));
pvals_sort = squeeze(pvals(o,:,:));

for i=1:782
ss(i,:) = minmax_norm(s(i,:));
ssc(i,:) = minmax_norm(sc(i,:));
end

for i=1:150
ss(:,i) = smooth(ss(:,i),40);
ssc(:,i) = smooth(ssc(:,i),40);
pv(:,i) = smooth(pvals_sort(:,i),40);
end
subplot(2,2,1)
imagesc(ss)
subplot(2,2,2)
imagesc(ssc)
subplot(2,2,3)
imagesc(pv)
