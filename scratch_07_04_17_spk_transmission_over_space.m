



% [ fields ] = bz_getPlaceFields1D(rateMap{4});
trialINts = behavior.events.trialIntervals(behavior.events.trialConditions==4,:);
for i=54:143
sp1 = Restrict(spikes.times{i},trialINts);
sp_ls = Restrict(spikes.times{14},trialINts);
sp_phase = Restrict([spikes.times{i} ,spkphases{i}'],trialINts);
sp_phase = sp_phase(:,2);
for j=1:length(bins)-1
    f = find(sp_phase>bins(j));
    ff = find(sp_phase<bins(j+1));
    fff=  intersect(f,ff);
    s{1} = sp_ls;
    s{2} = sp1(fff);
    if ~isempty(s{2})
        [times groups]=spikes2sorted(s);
        [ccg t] = CCG(times,groups,'binSize',.0005,'duration',.2);
        crosscorrs(i,j,:) = ccg(:,1,2);
    end
    end
end


for i=50:143
imagesc(t,bins,squeeze(crosscorrs(i,:,:)))
if ~isempty(fields{i})
   title(fields{i}{1}.COM)
end
pause
end


c=1;clear cc bb
for i=50:143
imagesc(t,bins,squeeze(crosscorrs(i,:,:)))
if ~isempty(fields{i})
title(fields{i}{1}.COM)
cc(c) = fields{i}{1}.COM;
bb(c) = i;
c=1+c;
end
end

[a b]=sort(cc);
for i=1:length(b)
subplot(4,3,i)
imagesc(t,bins,squeeze(crosscorrs(bb(b(i)),:,:)))
subplot(4,3,i)
title(cc(b(i)))
end