clear s p pos

for i=1:length(spk_trains{2})
s{i} = makeLength(spk_trains{2}{i}(2,:),200);
p{i} = makeLength(phase_trains{2}{i}(2,:),200);
pos{i} = makeLength(position{2},200);
end



ss = cell2mat(s');
pp = cell2mat(p');
ppp(1,:,:) = pp;
sss(1,:,:) = ss;

for smoothing = 14:1000
[track_info,pos_info_val] = Info_Analysis(discretize(sss,30),1,smoothing);
info(smoothing+1,:) = track_info;
[track_info,pos_info_val] = Info_Analysis(discretize(ppp,30),1,smoothing);
infop(smoothing+1,:) = track_info;
subplot(2,2,1)
imagesc(info)
subplot(2,2,2)
imagesc(infop)
subplot(2,2,3)
plot(smoothing,max(infop(smoothing,:)),'.k')
hold on
plot(smoothing,max(info(smoothing,:)),'.r')
subplot(2,2,4)  

pause(.01)
end


for i=1:10:3000
t = circ_smoothTS(p(2,:),i,'method','mean');
tt = smooth(s(2,:),i);
for iter = 1:5
r = randperm(length(s));
[b  dev]=glmfit([cos(t(r(1:length(r)./2))) sin(t(r(1:length(r)./2))) t(r(1:length(r)./2))],pos(r(1:length(r)./2)),'normal');
yfit = glmval(b,[cos(t(r(length(r)./2+1:end))) sin(t(r(length(r)./2+1:end))) t(r(length(r)./2+1:end))],'identity');
mse_phase(iter) = mean((yfit-pos(r(length(r)./2+1:end))').^2);
[b  dev]=glmfit([(tt(r(1:length(r)./2)))],pos(r(1:length(r)./2)),'normal');
yfit = glmval(b,tt(r(length(r)./2+1:end)),'identity');
mse_rate(iter) = mean((yfit-pos(r(length(r)./2+1:end))').^2);
end
plot(i,mean(mse_rate),'.r')
hold on
plot(i,mean(mse_phase),'.g')
pause(.01)
end
