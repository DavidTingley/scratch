d = dir('*201*');

c=1;
for i=1:length(d)
    cd(d(i).name)
    spikes = bz_GetSpikes('noprompts',true,'region','hpc');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    state = bz_LoadStates(pwd,'SleepState');
    if ~isempty(spikes) & ~isempty(ripples) & isfield(state.ints,'NREMstate')
        if length(ripples.peaks) > 15
        [times groups] = spikes2sorted([Restrict(ripples.peaks,double(state.ints.NREMstate)),spikes.times]);
        [ccg t] = CCG(times,groups,'binSize',.001);
        for j = 1:length(spikes.times)
        cc_hpc_nrem(c,:) = ccg(:,1,j);c=1+c;
        end
        end
    end
    cd /home/david/datasets/ripples_LS/
end
c=1
for i=1:length(d)
    cd(d(i).name)
    spikes = bz_GetSpikes('noprompts',true,'region','ls');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    state = bz_LoadStates(pwd,'SleepState');
    if ~isempty(spikes) & ~isempty(ripples) & isfield(state.ints,'NREMstate')
        if length(ripples.peaks) > 15
        [times groups] = spikes2sorted([Restrict(ripples.peaks,double(state.ints.NREMstate)),spikes.times]);
        [ccg t] = CCG(times,groups,'binSize',.001);
        for j = 1:length(spikes.times)
        cc_ls_nrem(c,:) = ccg(:,1,j);c=1+c;
        end
        end
    end
    cd /home/david/datasets/ripples_LS/
end

c=1;
for i=1:length(d)
    cd(d(i).name)
    spikes = bz_GetSpikes('noprompts',true,'region','hpc');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    state = bz_LoadStates(pwd,'SleepState');
    if ~isempty(spikes) & ~isempty(ripples)
        if length(ripples.peaks) > 15
        [times groups] = spikes2sorted([Restrict(ripples.peaks,double(state.ints.WAKEstate)),spikes.times]);
        [ccg t] = CCG(times,groups,'binSize',.001);
        for j = 1:length(spikes.times)
        cc_hpc_wake(c,:) = ccg(:,1,j);c=1+c;
        end
        end
    end
    cd /home/david/datasets/ripples_LS/
end
c=1
for i=1:length(d)
    cd(d(i).name)
    spikes = bz_GetSpikes('noprompts',true,'region','ls');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    state = bz_LoadStates(pwd,'SleepState');
    if ~isempty(spikes) & ~isempty(ripples)
        if length(ripples.peaks) > 15
        [times groups] = spikes2sorted([Restrict(ripples.peaks,double(state.ints.WAKEstate)),spikes.times]);
        [ccg t] = CCG(times,groups,'binSize',.001);
        for j = 1:length(spikes.times)
        cc_ls_wake(c,:) = ccg(:,1,j);c=1+c;
        end
        end
    end
    cd /home/david/datasets/ripples_LS/
end

c=1
for i=1:length(d)
    cd(d(i).name)
    spikes = bz_GetSpikes('noprompts',true,'region','ls');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    state = bz_LoadStates(pwd,'SleepState');
    if ~isempty(spikes) & ~isempty(ripples)
        if length(ripples.peaks) > 15
        [times groups] = spikes2sorted([ripples.peaks,spikes.times]);
        [ccg t] = CCG(times,groups,'binSize',.001);
        for j = 1:length(spikes.times)
        cc_ls(c,:) = ccg(:,1,j);c=1+c;
        end
        end
    end
    cd /home/david/datasets/ripples_LS/
end


whos cc*