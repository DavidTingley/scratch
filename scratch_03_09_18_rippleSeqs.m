
% d = dir('*201*');
% for rec=1:length(d)
%     cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    if ~isempty(ripples)
        spikes = bz_GetSpikes;
        for i=1:length(spikes.times)
            if strcmp(spikes.region{i},'ls')
                reg(i) = 1;
            elseif strcmp(spikes.region{i},'ca1') | strcmp(spikes.region{i},'hpc') | strcmp(spikes.region{i},'hpc')
                reg(i) = 2;
            else
                reg(i) = nan;
            end
        end
        numCells = length(spikes.times);

        %% runa actual data
        parfor i=1:length(ripples.timestamps)
            spk{i} = zeros(numCells,101);
            for ii=1:numCells
                f = Restrict(spikes.times{ii},[ripples.peaks(i)-.049999 ripples.peaks(i)+.049999]);
                spk{i}(ii,ceil(51+(f-ripples.peaks(i))*1000))=1;
                spk_smooth{i}(ii,:) = [zeros(25,1); Smooth(spk{i}(ii,:),5); zeros(25,1)];
            end
        end
    
        dat = cell2mat(spk);
        dat_smooth = cell2mat(spk_smooth);
        
        
        idx = find(reg==1); % LS
        [W_ls, H_ls, cost_ls,loadings_ls,power_ls] = seqNMF(dat_smooth(idx,:),'L',100,'K',40,'lambda',.000001,'showplot',0);
%         clear ind spk spk_smooth dat dat_smooth
        
        idx = find(reg==2); % HPC
        [W_hpc, H_hpc, cost_hpc,loadings_hpc,power_hpc] = seqNMF(dat_smooth(idx,:),'L',100,'K',40,'lambda',.000001,'showplot',0);
%         clear ind spk spk_smooth dat dat_smooth

%% run shuffled data..
        for iter = 1:20
         parfor i=1:length(ripples.timestamps)
            spk_shuffle{i} = zeros(numCells,101);
            r = randperm(length(ripples.peaks)); % random
            for ii=1:numCells
                f = Restrict(spikes.times{ii},[ripples.peaks(r(i))-.049999 ripples.peaks(r(i))+.049999]);
                spk_shuffle{i}(ii,ceil(51+(f-ripples.peaks(r(i)))*1000))=1;
                spk_smooth_shuffled{i}(ii,:) = [zeros(25,1); Smooth(spk_shuffle{i}(ii,:),5); zeros(25,1)];
            end
         end
        dat_smooth_shuffle = cell2mat(spk_smooth_shuffled);
        idx = find(reg==2); % HPC
        [W_hpc_shuffle, H_hpc_shuffle, cost_hpc_shuffle,loadings_hpc_shuffle(iter,:),power_hpc_shuffle(iter)] = seqNMF(dat_smooth_shuffle(idx,:),'L',100,'K',40,'lambda',.000001,'showplot',0);
        idx = find(reg==1); % LS
        [W_ls_shuffle, H_ls_shuffle, cost_ls_shuffle,loadings_ls_shuffle(iter,:),power_ls_shuffle(iter)] = seqNMF(dat_smooth_shuffle(idx,:),'L',100,'K',40,'lambda',.000001,'showplot',0);
        end
    end
    
    save([sessionInfo.FileName '.seqNMF.mat'],'-v7.3')
    
% end