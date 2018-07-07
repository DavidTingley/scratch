
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


        parfor i=1:length(ripples.timestamps)
            spk{i} = zeros(numCells,101);
            for ii=1:numCells
                f = Restrict(spikes.times{ii},[ripples.peaks(i)-.049999 ripples.peaks(i)+.049999]);
                spk{i}(ii,ceil(51+(f-ripples.peaks(i))*1000))=1;

                spk_smooth{i}(ii,:) = [zeros(25,1); Smooth(spk{i}(ii,:),5); zeros(25,1)];
%             if strcmp(spikes.region{ii},'ls')% | strcmp(spikes.region{ii},'ca1') | strcmp(spikes.region{ii},'ca3')
%                 ind(ii) = 1;
%             end
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
    end
    
    save([sessionInfo.FileName '.seqNMF.mat'],'-v7.3')
    
% end