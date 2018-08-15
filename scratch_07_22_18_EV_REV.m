% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');
binSize = [.003 .005 .01 .025 .05 .1 .25 .5 1 10];


for ii=1:95
    cd(d(ii).name)
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    if ~isempty(ls_spikes) & ~isempty(SleepState.ints.NREMstate) & exist([ls_spikes.sessionName '.behavior.mat'])
    if ~isempty(ls_spikes.times)
    ls_spikesNREM = ls_spikes;
    load([ls_spikes.sessionName '.behavior.mat'])
    for spk = 1:length(ls_spikes.times)
%        ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},double(SleepState.ints.NREMstate)); 
       ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},behavior.events.trialIntervals(behavior.events.trialConditions==2,:));
%        ls_spikesNREM.times{spk} = Restrict(ls_spikes.times{spk},behavior.events.trialIntervals); 
    end
    
   
    for bins = 1:length(binSize)
   
    spkmat_ls = bz_SpktToSpkmat(ls_spikes.times,'binSize',binSize(bins));
    spkmatNREM_ls = bz_SpktToSpkmat(ls_spikesNREM.times,'binSize',binSize(bins));

    for i=1:length(behavior.events.trials)
       bStart(i) = behavior.events.trials{i}.timestamps(1);
       bStop(i) = behavior.events.trials{i}.timestamps(end); 
    end
    pre = [0 min(bStart)];
    post = [max(bStop) spkmat_ls.timestamps(end)];
    behav = [pre(end) post(1)]; clear bStop bStart
    
    %% get initial correlations
    intervals = [pre; behav; post];
    % [corrmat pval] = corr(spkmat.data);

    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmat_ls.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmat_ls.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmat_ls.data(start:stop,:));
        if int==2 % int == 1 | int == 3
            [nah start] = min(abs(spkmatNREM_ls.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_ls.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_ls.data(start:stop,:));
        end
    end 

    for i = 1:size(intervals,1)
        for j=1:size(intervals,1)
            temp = corrcoef(corrmat{i},corrmat{j},'rows','complete');
            Rmod(i,j) = temp(2);
        end
    end

    EV_ls(bins,ii) = ((Rmod(2,3) - Rmod(2,1) * Rmod(1,3)) / sqrt((1 - Rmod(2,1)^2) * (1 - Rmod(1,3)^2)))^2;
    REV_ls(bins,ii) = ((Rmod(2,1) - Rmod(2,3) * Rmod(1,3)) / sqrt((1 - Rmod(2,3)^2) * (1 - Rmod(1,3)^2)))^2;


    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    if ~isempty(hpc_spikes.times)
    hpc_spikesNREM = hpc_spikes;
    for spk = 1:length(hpc_spikes.times)
%        hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},double(SleepState.ints.NREMstate));   
       hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals(behavior.events.trialConditions==2,:)); 
%        hpc_spikesNREM.times{spk} = Restrict(hpc_spikes.times{spk},behavior.events.trialIntervals); 
    end

    spkmat_hpc = bz_SpktToSpkmat(hpc_spikes.times,'binSize',binSize(bins));
    spkmatNREM_hpc = bz_SpktToSpkmat(hpc_spikesNREM.times,'binSize',binSize(bins));

    for i=1:length(behavior.events.trials)
       bStart(i) = behavior.events.trials{i}.timestamps(1);
       bStop(i) = behavior.events.trials{i}.timestamps(end); 
    end
    pre = [0 min(bStart)];
    post = [max(bStop) spkmat_hpc.timestamps(end)];
    behav = [pre(end) post(1)]; clear bStop bStart
    %% get initial correlations
    intervals = [pre; behav; post];
    % [corrmat pval] = corr(spkmat.data);

    for int = 1:size(intervals,1)
        [nah start] = min(abs(spkmat_hpc.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmat_hpc.timestamps-intervals(int,2)));

        [corrmat{int} pval{int}] = corr(spkmat_hpc.data(start:stop,:));
        if int==2 % int == 1 | int == 3
            [nah start] = min(abs(spkmatNREM_hpc.timestamps-intervals(int,1)));
            [nah stop] = min(abs(spkmatNREM_hpc.timestamps-intervals(int,2)));
            [corrmat{int} pval{int}] = corr(spkmatNREM_hpc.data(start:stop,:));
        end
    end

    for i = 1:size(intervals,1)
        for j=1:size(intervals,1)
            temp = corrcoef(corrmat{i},corrmat{j},'rows','complete');
            Rmod(i,j) = temp(2);
        end
    end

    EV_hpc(bins,ii) = ((Rmod(2,3) - Rmod(2,1) * Rmod(1,3)) / sqrt((1 - Rmod(2,1)^2) * (1 - Rmod(1,3)^2)))^2;
    REV_hpc(bins,ii) = ((Rmod(2,1) - Rmod(2,3) * Rmod(1,3)) / sqrt((1 - Rmod(2,3)^2) * (1 - Rmod(1,3)^2)))^2;
    end
    
    %% cross region
    spikes = bz_GetSpikes('noprompts',true);
    spikesNREM = spikes;
    for spk = 1:length(spikes.times)
%        spikesNREM.times{spk} = Restrict(spikes.times{spk},double(SleepState.ints.NREMstate));   
        spikesNREM.times{spk} = Restrict(spikes.times{spk},behavior.events.trialIntervals(behavior.events.trialConditions==2,:));  
%         spikesNREM.times{spk} = Restrict(spikes.times{spk},behavior.events.trialIntervals); 
       for spk2 = 1:length(spikes.times)
           if strcmp(spikes.region{spk},'hpc') & strcmp(spikes.region{spk2},'ls')
              idx(spk,spk2) = 1;
           else
               idx(spk,spk2) = 0;
           end
       end
    end

    spkmat = bz_SpktToSpkmat(spikes.times,'binSize',binSize(bins));
    spkmatNREM = bz_SpktToSpkmat(spikesNREM.times,'binSize',binSize(bins));
    
    for int = 1:size(intervals,1)
    [nah start] = min(abs(spkmat.timestamps-intervals(int,1)));
    [nah stop] = min(abs(spkmat.timestamps-intervals(int,2)));

    [corrmat{int} pval{int}] = corr(spkmat.data(start:stop,:));
    if int==2 % int == 1 | int == 3
        [nah start] = min(abs(spkmatNREM.timestamps-intervals(int,1)));
        [nah stop] = min(abs(spkmatNREM.timestamps-intervals(int,2)));
        [corrmat{int} pval{int}] = corr(spkmatNREM.data(start:stop,:));
    end
        corrmat{int}(find(idx==0)) = nan;
    end
    clear idx

    for i = 1:size(intervals,1)
        for j=1:size(intervals,1)
            temp = corrcoef(corrmat{i},corrmat{j},'rows','complete');
            Rmod(i,j) = temp(2);
        end
    end

    EV_cross(bins,ii) = ((Rmod(2,3) - Rmod(2,1) * Rmod(1,3)) / sqrt((1 - Rmod(2,1)^2) * (1 - Rmod(1,3)^2)))^2;
    REV_cross(bins,ii) = ((Rmod(2,1) - Rmod(2,3) * Rmod(1,3)) / sqrt((1 - Rmod(2,3)^2) * (1 - Rmod(1,3)^2)))^2;
    
    
    behavType{ii} = behavior.description;
    EV_ls(EV_ls==0) = nan;
    EV_hpc(EV_hpc==0) = nan;
    REV_ls(REV_ls==0) = nan;
    REV_hpc(REV_hpc==0) = nan;
    EV_cross(EV_cross==0) = nan;
    REV_cross(REV_cross==0) = nan;
    
    subplot(3,4,bins)
    histogram(EV_ls(bins,EV_ls(bins,:)~=0)-REV_ls(bins,EV_ls(bins,:)~=0),-1:.1:1,'normalization','pdf','FaceColor','m')
    hold on
    histogram(EV_hpc(bins,EV_hpc(bins,:)~=0)-REV_hpc(bins,EV_hpc(bins,:)~=0),-1:.1:1,'normalization','pdf','FaceColor','k')
    hold off
    subplot(3,4,10)
    plot(nanmedian(EV_cross,2),'.b')
    hold on
    plot(nanmedian(REV_cross,2),'.r')
    hold off
    subplot(3,4,11)
    plot(nanmedian(EV_ls,2),'.b')
    hold on
    plot(nanmedian(REV_ls,2),'.r')
    title('ls')
    hold off
    subplot(3,4,12)
    plot(nanmedian(EV_hpc,2),'.b')
    hold on
    plot(nanmedian(REV_hpc,2),'.r')
    hold off
    title('hpc')

    pause(.1)

    end
    end
    end
cd ~/datasets/ripples_LS/
end