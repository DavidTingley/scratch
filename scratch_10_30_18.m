d = dir('*201*');
c=1;

for i=1:length(d)
    cd(d(i).name)
    
    rips = bz_LoadEvents(pwd,'LSRipples');
    
    if ~isempty(rips)
       ls_spikes = bz_GetSpikes('noprompts',true,'region','ls');
       
       if ~isempty(ls_spikes)
           sessionInfo = bz_getSessionInfo;

           for spk = 1:length(ls_spikes.times)
              ls_spikes.times{spk} = Restrict(ls_spikes.times{spk},[rips.timestamps(:,1)-.025 rips.timestamps(:,2)+.025]); 
           end
           [times groups] = spikes2sorted(ls_spikes.times);

           [ccg t] = CCG(times,groups,'duration',.005,'binSize',.000025); %.001
            offsets = ls_spikes.shankID-ls_spikes.shankID';


           for spk = 1:size(ccg,2)
               for spk2 = spk:size(ccg,3)
                   if spk ~= spk2
                   ccgs(c,:) = ccg(:,spk,spk2);
                   offs(c) = offsets(spk,spk2);
                   animal(c) = sum(double(sessionInfo.animal));
                   c=1+c;
                   end
               end
           end
       end
    end    
    cd /home/david/datasets/ripples_LS 
end

u = unique(animal);
asym = sum(ccgs(:,47:51)')-sum(ccgs(:,51:55)');
for i=1:size(ccgs,1)
    ccgs_z(i,:) = zscore(ccgs(i,:));
end

for i=1:length(u)
    subplot(3,2,i)    
    for o=-8:8
%           [pk loc] = max(nanmean(ccgs_z(offs==o & u(i) == animal,:)));
%           plot(o,loc,'.k')
%           hold on
        errorbar(o,nanmean(asym(offs==o & u(i) == animal)),sem(asym(offs==o & u(i) == animal)'),'k')
        hold on
%           plot(zscore(nanmean(ccgs_z(offs==o & u(i) == animal,:))))
%           hold on
    end
    title(u(i))
end
