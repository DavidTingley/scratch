d = dir('*201*');
warning off

for threshold = 1%1:10
% c=1;
% for i=1:length(d)
%     cd(d(i).name)
%     spikes = bz_GetSpikes('noprompts',true,'region','hpc');
%     popBursts = bz_LoadEvents(pwd,'CA1popBursts');
%     state = bz_LoadStates(pwd,'SleepState');
%     if ~isempty(spikes) & ~isempty(popBursts) & isfield(state.ints,'NREMstate')
%         popBursts = bz_FindPopBursts(spikes,'threshold',threshold,'durations',[10 500]);
%         popBursts.bursts = popBursts.bursts;
%         if length(popBursts.bursts) > 15
%         [times groups] = spikes2sorted([Restrict(popBursts.bursts,double(state.ints.NREMstate)),spikes.times]);
%         [ccg t] = CCG(times,groups,'binSize',.01);
%         for j = 1:length(spikes.times)
%         cc_hpc_nrem(threshold,c,:) = ccg(:,1,j);c=1+c;
%         end
%         end
%     end
% %     cd /home/david/datasets/popBursts_LS/
% cd E:\datasets\popBursts_LS
% end
c = 1;
cc = 1;
for i=1:length(d)
    cd(d(i).name)
    spikes = bz_GetSpikes('noprompts',true,'region','ls');
    hpc = bz_GetSpikes('noprompts',true,'region','hpc');
    sessionInfo =bz_getSessionInfo;
    state = bz_LoadStates(pwd,'SleepState');
    if ~isempty(spikes)  & isfield(state.ints,'NREMstate') & exist([sessionInfo.FileName '.CA1Ripples.events.mat'])
%         popBursts = bz_FindPopBursts(hpc,'threshold',threshold,'durations',[10 500],'binSize',.01);
        popBursts = bz_LoadEvents(pwd,'popBursts');
        ripples = bz_LoadEvents(pwd,'CA1Ripples');
        ripples_ls = bz_LoadEvents(pwd,'LSRipples');
        popBursts.bursts = ripples.peaks;
        for evt = 1:length(ripples.peaks)
           nSpikes{i}(evt) = length(Restrict(spikes.spindices(:,1),[ripples.peaks(evt)-.025 ripples.peaks(evt)+.025]));             
        end
        amplitudes{i} = zscore(ripples.data.peakAmplitude);
        durations{i} = zscore(ripples.data.duration);
        nCells(i) = length(spikes.times);
        
        if ~isempty(popBursts) & ~isempty((state.ints.NREMstate))
        %% NREM
%         [times groups] = spikes2sorted([Restrict(popBursts.bursts,double(state.ints.NREMstate)),spikes.times]);
%         [ccg t] = CCG(times,groups,'duration',1,'binSize',.001,'norm','rate');
%         for j = 1:length(spikes.times)
%             cc_ls_nrem(threshold,c,:) = (ccg(:,1,j));
%             c=1+c;
%         end
        
%         ls_bursts = bz_FindPopBursts(spikes,'threshold',3,'durations',[10 500],'binSize',.01);
        ls_bursts.bursts = ripples_ls.peaks;
        if ~isempty(ls_bursts.bursts)
        [times groups] = spikes2sorted({Restrict(popBursts.bursts,double(state.ints.NREMstate)),Restrict(ls_bursts.bursts,double(state.ints.NREMstate))});
        if ~isempty(times)
            [ccg t] = CCG(times,groups,'binSize',.001);
            cc_ls_nrem_pop(threshold,i,:) = (ccg(:,1,2));
        end
        end
        
        
        %% WAKE
%         [times groups] = spikes2sorted([Restrict(popBursts.bursts,double(state.ints.WAKEstate)),spikes.times]);
%         [ccg t] = CCG(times,groups,'duration',1,'binSize',.001,'norm','rate');
%         for j = 1:length(spikes.times)
%             cc_ls_wake(threshold,cc,:) = (ccg(:,1,j));
%             cc=1+cc;
%         end
        if ~isempty(ls_bursts.bursts)
        [times groups] = spikes2sorted({Restrict(popBursts.bursts,double(state.ints.WAKEstate)),Restrict(ls_bursts.bursts,double(state.ints.WAKEstate))});
        if ~isempty(times)
            [ccg t] = CCG(times,groups,'binSize',.001);
            cc_ls_wake_pop(threshold,i,:) = (ccg(:,1,2));
        end
        end
        end
    subplot(2,2,1)
    imagesc(squeeze(mean(cc_ls_nrem_pop,2)))
    caxis([-1 1])
    title('LS nrem popBursts')
    subplot(2,2,2)
    imagesc(squeeze(mean(cc_ls_wake_pop,2)))
    caxis([-1 1])
    title('LS wake popBursts')
%     subplot(2,2,4)
%     imagesc(squeeze(mean(cc_ls_wake,2)))
%     caxis([-2 2])
%     title('LS wake CCGs')
%     subplot(2,2,3)
%     imagesc(squeeze(mean(cc_ls_nrem,2)))
%     caxis([-2 2])
%     title('LS nrem CCGs')
    pause(.01)
    end
    cd /home/david/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
end

% c=1;
% for i=1:length(d)
%     cd(d(i).name)
%     spikes = bz_GetSpikes('noprompts',true,'region','hpc');
%     popBursts = bz_LoadEvents(pwd,'CA1popBursts');
%     state = bz_LoadStates(pwd,'SleepState');
%     if ~isempty(spikes) & ~isempty(popBursts)
%         popBursts = bz_FindPopBursts(spikes,'threshold',threshold,'durations',[10 500]);
%         popBursts.bursts = popBursts.bursts;
%         if length(popBursts.bursts) > 15
%         [times, groups] = spikes2sorted([Restrict(popBursts.bursts,double(state.ints.WAKEstate)),spikes.times]);
%         [ccg t] = CCG(times,groups,'binSize',.01);
%         for j = 1:length(spikes.times)
%         cc_hpc_wake(threshold,c,:) = ccg(:,1,j);c=1+c;
%         end
%         end
%     end
% %     cd /home/david/datasets/ripples_LS/
%     cd E:\datasets\ripples_LS
% end
% cc=1
% for i=1:length(d)
%     cd(d(i).name)
%     spikes = bz_GetSpikes('noprompts',true,'region','ls');
%     hpc = bz_GetSpikes('noprompts',true,'region','hpc');
%     popBursts = bz_LoadEvents(pwd,'CA1popBursts');
%     state = bz_LoadStates(pwd,'SleepState');
%     if ~isempty(spikes) & ~isempty(popBursts) & ~isempty(hpc)
%         popBursts = bz_FindPopBursts(hpc,'threshold',threshold,'durations',[10 500]);
%         popBursts.bursts = popBursts.bursts;
%         if length(popBursts.bursts) > 15
%         [times groups] = spikes2sorted([Restrict(popBursts.bursts,double(state.ints.WAKEstate)),spikes.times]);
%         [ccg t] = CCG(times,groups,'binSize',.01);
%         for j = 1:length(spikes.times)
%         cc_ls_wake(threshold,cc,:) = ccg(:,1,j);cc=1+cc;
%         end
%         
%         ls_bursts = bz_FindPopBursts(spikes,'threshold',3,'durations',[10 500]);
%         [times groups] = spikes2sorted({Restrict(popBursts.bursts,double(state.ints.WAKEstate)),Restrict(ls_bursts.bursts,double(state.ints.WAKEstate))});
%         if ~isempty(times)
%             [ccg t] = CCG(times,groups,'binSize',.01);
%             cc_ls_wake_pop(threshold,i,:) = zscore(ccg(:,1,2));
%         end
%         end
%     subplot(2,2,2)
%     imagesc(squeeze(mean(cc_ls_wake_pop,2)))
%     title('LS wake')
%     caxis([-1 1])
%     pause(.01)
%     end
%     cd /home/david/datasets/ripples_LS/
% % cd E:\datasets\ripples_LS
% end

% c=1
% for i=1:length(d)
%     cd(d(i).name)
%     spikes = bz_GetSpikes('noprompts',true,'region','ls');
%     popBursts = bz_LoadEvents(pwd,'CA1popBursts');
%     state = bz_LoadStates(pwd,'SleepState');
%     if ~isempty(spikes) & ~isempty(popBursts)
%         popBursts = bz_FindPopBursts(spikes,'threshold',threshold,'durations',[10 500]);
%         popBursts.bursts = popBursts.bursts;
%         if length(popBursts.bursts) > 15
%         [times groups] = spikes2sorted([popBursts.bursts,spikes.times]);
%         [ccg t] = CCG(times,groups,'binSize',.01);
%         for j = 1:length(spikes.times)
%         cc_ls(threshold,c,:) = ccg(:,1,j);c=1+c;
%         end
%         end
%     end
% %     cd /home/david/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
% end
% 
% c=1
% for i=1:length(d)
%     cd(d(i).name)
%     spikes = bz_GetSpikes('noprompts',true,'region','hpc');
%     popBursts = bz_LoadEvents(pwd,'CA1popBursts');
%     state = bz_LoadStates(pwd,'SleepState');
%     if ~isempty(spikes) & ~isempty(popBursts)
%         popBursts = bz_FindPopBursts(spikes,'threshold',threshold,'durations',[10 500]);
%         popBursts.bursts = popBursts.bursts;
%         if length(popBursts.bursts) > 15
%         [times groups] = spikes2sorted([popBursts.bursts,spikes.times]);
%         [ccg t] = CCG(times,groups,'binSize',.01);
%         for j = 1:length(spikes.times)
%         cc_hpc(threshold,c,:) = ccg(:,1,j);c=1+c;
%         end
%         end
%     end
% %     cd /home/david/datasets/ripples_LS/
%     cd E:\datasets\ripples_LS
% end

end
whos cc*