


for i=1:140
    sp{i}=[];
    sp_all{i} = [];
    for j=1:length(behavior.events.trialIntervals)
        if behavior.events.trialConditions(j)==2
        if ~strcmp(spikes.region{i},'ls')
        sp{i} = [sp{i}; Restrict(spikes.times{i},behavior.events.trialIntervals(j,:))];
        end
        sp_all{i} = [sp_all{i}; Restrict(spikes.times{i},behavior.events.trialIntervals(j,:))];
        end
    end
end

[tree,counts,siz] = EventTree(sp,.035,2,8);
f = find(tree{end}>10);
clear ts

for i=1:length(f)
[seq] = sparse2mat(siz,f(i));
[t num{i}] = findchain(sp,seq,.035,8);

ts{i} = mean(t{1}')';
end

ts{end+1}=spikes.times{14};
[spk groups] = spikes2sorted(ts);
[cct t] = CCG(spk,groups,'binSize',.001);
[spk groups] = spikes2sorted(sp_all);
[ccg_all t] = CCG(spk,groups,'binSize',.001);


for i=1:length(f)
subplot(2,2,1)
plot(t,cct(:,i,end))
subplot(2,2,2)
plot(t,smooth(cct(:,i,end),25))
[seq] = sparse2mat(siz,f(i));
subplot(2,2,3)
for k=1:2
    plot(t,ccg_all(:,seq(k),14))
    hold on
end
hold off
subplot(2,2,4)
for k=1:2
    plot(t,zscore(smooth(ccg_all(:,seq(k),14),25)))
    hold on
end
hold off

title(seq)
pause
end