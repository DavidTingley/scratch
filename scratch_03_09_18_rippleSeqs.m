
% d = dir('*201*');
% for rec=1:length(d)
%     cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    if ~isempty(ripples)
        spikes = bz_GetSpikes;
        numCells = length(spikes.times);


        for i=1:length(ripples.timestamps)
            spk{i} = zeros(numCells,101);
            for ii=1:numCells
                f = Restrict(spikes.times{ii},[ripples.peaks(i)-.049999 ripples.peaks(i)+.049999]);
                spk{i}(ii,ceil(51+(f-ripples.peaks(i))*1000))=1;

                spk_smooth{i}(ii,:) = Smooth(spk{i}(ii,:),12);
            if strcmp(spikes.region{ii},'ls')% | strcmp(spikes.region{ii},'ca1') | strcmp(spikes.region{ii},'ca3')
                ind(ii) = 1;
            end
            end
        end
    
        dat = cell2mat(spk);
        dat_smooth = cell2mat(spk_smooth);
        ind = find(ind);
        
        [W, H, cost,loadings,power] = seqNMF(dat_smooth(ind,:),'L',100,'K',40,'lambda',.1,'showplot',0);
        clear ind spk spk_smooth dat dat_smooth
    end
    
% end