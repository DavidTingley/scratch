function [] = scratch_07_22_18_bayesianTemplate()
% cd D:\Dropbox\datasets\lsDataset
d = dir('*201*');
binSize = [.01];
count_ls = 1;
count_hpc = 1;
overlap = 1;


% load('D:\Dropbox\Documents\pubs\inProgress\2019_ReplayDavid\data\20170429_0um_0um_merge/exampleOrder.mat', 'b','ex')
% seqOrder = b;

    sessionInfo = bz_getSessionInfo;
    disp(['working on ' sessionInfo.FileName])
    ls_spikes = bz_GetSpikes('region','ls','noprompts',true);
    hpc_spikes = bz_GetSpikes('region','hpc','noprompts',true);
    spikes = bz_GetSpikes('noprompts',true);
    SleepState = bz_LoadStates(pwd,'SleepState');
    ripples = bz_LoadEvents(pwd,'CA1Ripples');

    
    if exist([sessionInfo.FileName '.firingMaps.cellinfo.mat']) 
    
    ls_spikesNREM = ls_spikes;
    ls_spikesBEHAV = ls_spikes;
    idx_ls = find(strcmp(spikes.region,'ls'));
    
%     hpcRegion{count_hpc} = 'ca1';
%     if isempty(hpc_spikes) % comment out to exclude CA3 recs
%         hpc_spikes = bz_GetSpikes('region','ca3','noprompts',true);
%         hpcRegion{count_hpc} = 'ca3';
%     end
    idx_hpc = find(strcmp(spikes.region,'hpc') | strcmp(spikes.region,'ca3'));
    
    if exist([sessionInfo.FileName '.behavior.mat']) & ...
        ~isempty(ripples) & ...
        isfield(SleepState.ints,'NREMstate') & ~isempty(SleepState.ints.NREMstate)
        load([sessionInfo.FileName '.behavior.mat'])
        load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
        [firingMaps] = bz_firingMap1D(spikes,behavior,4);
    
    %%
%     lfp = bz_GetLFP([95]);
    lfp = bz_GetLFP([sessionInfo.ca1]);
    seqOrder = 1:length(spikes.times);
    ex = zeros(length(spikes.times),1);
    %%  

    for bins = 1:length(binSize)
%         hpc_spikesNREM = hpc_spikes;
%         for spk = 1:length(spikes.times) % restricting speeds up spkMat processing
%            hpc_spikesNREM.times{spk} = Restrict(spikes.times{spk},[ripples.peaks-.5 ripples.peaks+.5]); 
%         end
        
        spkmatNREM_hpc = bz_SpktToSpkmat(spikes.times,'overlap',overlap,'binSize',binSize(bins)* overlap);
        nCells = nan(length(firingMaps.rateMaps),length(ripples.peaks));
        nSpks = nan(length(firingMaps.rateMaps),length(ripples.peaks));
        rankOrder = nan(length(firingMaps.rateMaps),length(ripples.peaks));
        outR = nan(length(firingMaps.rateMaps),length(ripples.peaks));
        outR_shuf = nan(length(firingMaps.rateMaps),length(ripples.peaks),100);
        integral_hpc = nan(length(firingMaps.rateMaps),length(ripples.peaks));
        integral_hpc_shuf = nan(length(firingMaps.rateMaps),length(ripples.peaks),100);
        corrs_hpc = nan(length(firingMaps.rateMaps),length(ripples.peaks));
        corrs_hpc_shuf = nan(length(firingMaps.rateMaps),length(ripples.peaks),100);
        slope_hpc =  nan(length(firingMaps.rateMaps),length(ripples.peaks));
        slope_hpc_shuf =  nan(length(firingMaps.rateMaps),length(ripples.peaks));

        for t = 1:length(firingMaps.rateMaps)
            if size(firingMaps.rateMaps{t},2) >= 9  % require minimum of 9 good trials to look for replay of it

        
             for spk=1:length(spikes.times)
                pf(spk) = ~isempty(fields{t}{spk}); %% all place fields
%                 pf(spk) = length(fields{t}{spk})==1;  %% only one place field per cell
             end
%              pf = find(pf);
             keep = intersect(find(pf),idx_hpc);%-length(idx_ls); % are you a place cell, in HPC?
%              keep = 1:length(idx_hpc);  % uncomment to keep all HPC cells that fired
%              keep(sum(rates')==0) = [];

            template = squeeze(mean(firingMaps.rateMaps{t},2));
            template = template(:,16:185);
%             for i = 1:size(template,1)
%                template(i,:) = mean_norm(Smooth(template(i,:),5)')';
%             end
            for event = 1:length(ripples.peaks)
                % start bigger than the event itself
%                 start = round((round(ripples.timestamps(event,1) * 1000)-50) ./ (spkmatNREM_hpc.dt*1000));
%                 stop = round((round(ripples.timestamps(event,2) * 1000)+50) ./ (spkmatNREM_hpc.dt*1000));
                start = round((round(ripples.peaks(event) * 1000)-50) ./ (spkmatNREM_hpc.dt*1000));
                stop = round((round(ripples.peaks(event) * 1000)+50) ./ (spkmatNREM_hpc.dt*1000));
                
                if  stop < size(spkmatNREM_hpc.data,1) & stop-start > 5
                    for spk = 1:size(spkmatNREM_hpc.data,2)
                            data(:,spk) = (spkmatNREM_hpc.data(start:stop,spk)')';  
                            counts(:,spk) = (spkmatNREM_hpc.data(start:stop,spk)')';  
                    end

                    % cut 0 rate bins..
                    while sum(counts(1,keep)) < 1 & size(counts,1) > 1
                       data = data(2:end,:);
                       counts = counts(2:end,:);
                    end
                    while sum(counts(end,keep)) < 1 & size(counts,1) > 1
                       data = data(1:end-1,:);
                       counts = counts(1:end-1,:);
                    end
                    
                    %% get max rate bin
%                     [a b] = max(sum(data,2));
    %                 or search from middle...
    %                 b = round(size(data,1)/2);
                     % find clipping point at beginning
    %                  sta = 1;
    %                 while sum(counts(b-sta,keep),2) > 1 & b-sta > 1 % max length is +/-50
    %                    sta = sta + 1;
    %                 end
    %                  % find clipping point at beginning
    %                  sto = 1;
    %                 while sum(counts(b+sto,keep),2) > 1 & sto +b < size(data,1)-1  % max length 500 ms
    %                    sto = sto + 1;
    %                 end

                    % redefine start/stop by pop burst..
    %                 data = data(b-sta:b+sto,:);


                if stop < size(spkmatNREM_hpc.data,1) & size(data,1) > 4 & sum(sum(counts(:,keep))>0) > 4
                    data = data ./ binSize;
                    [Pr, prMax] = placeBayes(data(:,keep), template(keep,:), spkmatNREM_hpc.dt);
                    
                    % now calc correlations and control distro w/ cut data...
                    [corrs_hpc(t,event) ]= corr([1:length(prMax)]',prMax,'rows','complete');
                    [outR(t,event) outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1));
                    nCells(t,event) = sum(sum(counts(:,keep))>0);
                    nSpks(t,event) = sum(sum(counts(:,keep)));

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    idx = intersect(find(nanmean(data)~=0),keep); % only take cells that spiked...               
                    [a b ord_template] = sort_cells(template(idx,:));
                    [a b ord_data] = sort_cells(data(:,idx)');
                    for i =1:size(data,2)
                        [ts(i)]=mean(find(data(:,i)>0));
                    end
                    [a b ord_avg] = sort_cells(ts(idx)'); clear ts
                    [rankOrder(t,event) pvals(t,event)] = corr(ord_template,ord_avg,'rows','complete');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                    for iter = 1:10
    %                     template_shuf = bz_shuffleCircular(squeeze(mean(firingMaps.rateMaps_unsmooth{t}(idx_hpc,:,:),2)));
    %                     for i = 1:size(template,1)
    %                        template_shuf(i,:) = mean_norm(Smooth(template_shuf(i,:),5)')';
    %                     end
                        template_shuf = bz_shuffleCellID(template);
                        [Pr_shuf prMax_shuf] = placeBayes((data(:,keep)')', template_shuf(keep,:), spkmatNREM_hpc.dt*150);
                        [outR_shuf(t,event,iter) outID] = makeBayesWeightedCorrBatch1(Pr_shuf,ones(size(Pr_shuf,1),1));
                        corrs_hpc_shuf(t,event,iter) = corr([1:length(prMax_shuf)]',prMax_shuf,'rows','complete');
                        
                        [a b ord_shuf] = sort_cells(template_shuf(idx,:));
                        [rankOrder_shuffle(t,event,iter)  pvals_shuf(t,event,iter)] = corr(ord_template,ord_shuf,'rows','complete');
                                              
                        [slope_hpc_shuf(t,event,iter) integral_hpc_shuf(t,event,iter)] = Pr2Radon(Pr_shuf');
                    end

                    if sum(sum(spkmatNREM_hpc.data(start:stop,:)))> 5 * overlap & sum(~isnan(sum(Pr')))>5
                        [slope_hpc(t,event) integral_hpc(t,event) ] = Pr2Radon(Pr');
                    else 
                        corrs_hpc(t,event) = nan;
                        corrs_hpc_shuf(t,event,1:100) = nan;
                        rankOrder(t,event) = nan;
                        outR(t,event) = nan;
                        outR_shuf(t,event,1:100) = nan;
                        slope_hpc(t,event) = nan;
                        integral_hpc(t,event) =nan;
                        slope_hpc_shuf(t,event) = nan;
                        integral_hpc_shuf(t,event,1:100) = nan;
                    end
                else
                    corrs_hpc(t,event) = nan;
                    corrs_hpc_shuf(t,event,1:100) = nan;
                    rankOrder(t,event) = nan;
                    outR(t,event) = nan;
                    outR_shuf(t,event,1:100) = nan;
                    slope_hpc(t,event) = nan;
                    integral_hpc(t,event) =nan;
                    slope_hpc_shuf(t,event) = nan;
                    integral_hpc_shuf(t,event,1:100) = nan;
                    Pr = [];
                end

                nCells(event) = length(keep);
                spkCount(event) = sum(sum(data(:,keep)))./overlap;
                eventDuration(event) = (stop-start)*spkmatNREM_hpc.dt;

    %% plotting
    if abs(outR(t,event))>.2
        interval = [ripples.timestamps(event,1)-.05 ripples.timestamps(event,2)+.05];
        subplot(4,2,1)
        bz_plotRasterTrial(spikes,interval,10,seqOrder,zeros(length(spikes.times),1));


        ylim([0 60])
        title(['R = ' num2str(corrs_hpc(t,event)) ', evt #' num2str(event)]);
        subplot(4,2,2)
        bz_plotEphys(lfp,'timewin',interval)


        rZscore = (abs(outR)-mean(abs(outR_shuf),3))./std(abs(outR_shuf),[],3);  % Grosmark/Buzsaki "sequence score"...
        subplot(4,2,3)
%         bz_plotRasterTrial(spikes,interval,10,seqOrder,((~ismember(spikes.UID,hpc_spikes.UID) + ex' + ~pf)==3)); 
        ylim([0 60])
        title(['condition: ' num2str(t) ', event: ' num2str(event)])
        hold off

        subplot(4,2,4) 
        % line([0 3],[.012 .012],'color','r');
        bz_plotEphys([],'spikes',spikes,...
                            'scalelfp',5,...
                            'plotcells',spikes.UID(keep),...
                            'timewin',[ripples.timestamps(event,1)-.05 ripples.timestamps(event,2)+.05])

        subplot(4,2,5)
        scatter(max((rZscore(:,1:event))),absmax(outR(:,1:event)),'.k')
        % scatter(absmax((corrs_hpc(:,1:event))),spkCount(1:event),'.k')
        xlabel('rZ radon integral')
        ylabel('weighted corr')
        title(['weighted corr:' num2str(outR(t,event))])

        subplot(4,2,6)
        title([num2str(integral_hpc(t,event)) ' / ' num2str(rZscore(t,event))])
        imagesc((Pr'))
        % 
        % title(corr(PR{count_hpc}(1:event)',max((d(:,1:event)))','rows','complete'))

        subplot(4,2,7)
        % plot(PR{count_hpc}(1:event),'.k')
        plot(absmax((corrs_hpc(:,1:event))),absmax(rZscore(:,1:event)),'.k')
        xlabel('rank order corr')
        ylabel('rZ score')
        % imagesc(data(:,ord)')
        title('rank ord corr vs linear weight corr')

        subplot(4,2,8)
        title(['radon int: ' num2str(nanmean(integral_hpc_shuf(t,event,:),3)) ' / ' num2str(rZscore(t,event))])
        imagesc(flipud(data(:,keep)'))
        
        
        pause(.01)
        clf
    end
            end
    clear data counts
            end
            end
            clear pf
        end


         behavType{count_hpc} = behavior.description;



        save([sessionInfo.FileName '.replayComparisonData_allFields.mat'],'-v7.3')
        end
    end

end
