clear s p pos

for i=1:21
s{i} = makeLength(spk_trains{5}{i}(80,:),3700);
p{i} = makeLength(phase_trains{5}{i}(80,:),3700);
pos{i} = makeLength(position{5},3700);
end



ss = cell2mat(s');
pp = cell2mat(p');
ppp(1,:,:) = pp;
sss(1,:,:) = ss;

for smoothing = 1:1000
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


