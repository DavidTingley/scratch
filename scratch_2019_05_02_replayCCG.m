% clear all
cd /home/david/datasets/ripples_LS
% cd /home/david/Dropbox/data
d = dir('*201*');
pre_excess_corr = [];
post_excess_corr = [];
pre_excess_jitter_corr = [];
post_excess_jitter_corr = [];
beh_excess_corr = [];
beh_excess_jitter_corr = [];

fieldCounts = [];
for w=1:250
    compression(w) = 240/(length(w:501-w+1));
end

for ii=1:98
    excess_pre{ii}=[];
    excess_post{ii}=[];
    excess_pre_jitter{ii}=[];
    excess_post_jitter{ii}=[];
    excess_max_pre{ii}=[];
    excess_max_post{ii}=[];
    excess_beh{ii}=[];
    excess_beh_jitter{ii}=[];
    
    cd(d(ii).name)


    spikes = bz_GetSpikes('noprompts',true);
    sessionInfo = bz_getSessionInfo;
    if ~isempty(spikes)
    ls_idx = strcmp(spikes.region,'ls');
%     ls_idx = strcmp(spikes.region,'hpc') | strcmp(spikes.region,'ca3');
    if ~isempty(ls_idx) & exist([sessionInfo.FileName '.LSRipples.events.mat']) & exist([sessionInfo.FileName '.behavior.mat'])
    load([sessionInfo.FileName '.LSRipples.events.mat'])
%     load([sessionInfo.FileName '.CA1Ripples.events.mat'])
    load([sessionInfo.FileName '.behavior.mat'])
    if exist([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
        load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
        nFields = zeros(length(spikes.times),1);
        for c=1:length(fields)
            for cc = 1:length(fields{c})
                if ~isempty(fields{c}{cc})
                    nFields(cc) = nFields(cc) + 1;
                end
            end
    end
    else
         nFields = zeros(length(spikes.times),1);
    end
    for i=1:length(spikes.times)
        for j = 1:length(spikes.times)
           sh(i,j) = spikes.shankID(i) == spikes.shankID(j);
        end
    spk.times{i} = Restrict(spikes.times{i},behavior.events.trialIntervals);
    end
    [times groups] = spikes2sorted(spk.times); clear spk
    [ccg_beh t] = CCG(times,groups,'binSize',.001,'duration',.5,'norm','rate');
    ccg_beh(end,length(spikes.times),length(spikes.times))=0;
    
    times_jitter = times + rand(length(times),1)./5;
    [ccg_beh_jitter t] = CCG(times_jitter,groups,'binSize',.001,'duration',.5,'norm','rate');
    ccg_beh_jitter(end,length(spikes.times),length(spikes.times))=0;
        
        
    ccg_behav{ii} = single(ccg_beh(:,ls_idx,ls_idx));
%     temporalBias_beh{ii} = nan(size(ccg_beh,2),size(ccg_beh,3));
    if ~isempty(times) & sum(sum(~sh(ls_idx,ls_idx))) > 0
        for i = 1:size(ccg_beh,2)
                for j = 1:size(ccg_beh,3)
                    if sum(ccg_beh(:,i,j)) > 5 & sh(i,j) == 0 & ls_idx(i) == 1 & ls_idx(j) == 1
    %                 b(i,j) = (sum(ccg_beh(1:100,i,j)))-(sum(ccg_beh(102:end,i,j)));
    %                 br(i,j) = (sum(ccg_beh(1:100,i,j)))-(sum(ccg_beh(102:end,i,j)));


                        temporalBias_beh{ii}(i,j) = single(sum(ccg_beh(51:250,i,j) - ccg_beh(252:451,i,j))) ./single(sum(ccg_beh(51:250,i,j) + ccg_beh(252:451,i,j))); 

                        sma = Smooth(ccg_beh(:,i,j),5);
                        sma = sma ./ sum(sma);
                        la = Smooth(ccg_beh(:,i,j),200);
                        la = la ./ sum(la);
                        excess_beh{ii}(i,j) = sma(251)-la(251);
                        excess_max_beh{ii}(i,j) = max(sma(225:275))-max(la(225:275));
                        
                        % jitter vers
                        sma = Smooth(ccg_beh_jitter(:,i,j),5);
                        sma = sma ./ sum(sma);
                        la = Smooth(ccg_beh_jitter(:,i,j),200);
                        la = la ./ sum(la);
                        excess_beh_jitter{ii}(i,j) = sma(251)-la(251);
                
                        nSpikes_beh{ii}(i,j) = sum(ccg_beh(:,i,j));
                        fieldComp{ii}(i,j) = double(nFields(i)>0) + double(nFields(j)>0);
                    else
                        fieldComp{ii}(i,j) = double(nFields(i)>0) + double(nFields(j)>0);
                        nSpikes_beh{ii}(i,j) = 0;
                        temporalBias_beh{ii}(i,j) = NaN;
                        excess_beh{ii}(i,j) = NaN;
                        excess_beh_jitter{ii}(i,j) = NaN;
                        excess_max_beh{ii}(i,j) = NaN;
                    end
                    
                    
                end
        end
    end
    idx = find(ripples.peaks>behavior.events.trialIntervals(end,1));
    if ~isempty(idx)
    for i=1:length(spikes.times)
    spk.times{i} = [Restrict(spikes.times{i},[ripples.timestamps(idx,1)-.05 ripples.timestamps(idx,2)+.05])];
%     spk.times{i} = [Restrict(spikes.times{i},[ripples.timestamps(idx(1),1)-.05 spikes.spindices(end,1)])];
    end
    [times groups] = spikes2sorted(spk.times); clear spk
    if ~isempty(times) & sum(sum(~sh(ls_idx,ls_idx))) > 0
        [ccg_ls_rip t] = CCG(times,groups,'binSize',.001,'duration',.5,'norm','rate');
        ccg_ls_rip(:,length(spikes.times),length(spikes.times))=0;
        
        times_jitter = times + rand(length(times),1)./5;
        [ccg_ls_rip_jitter t] = CCG(times_jitter,groups,'binSize',.001,'duration',.5,'norm','rate');
        ccg_ls_rip_jitter(end,length(spikes.times),length(spikes.times))=0;
        
        ccg_post{ii} = single(ccg_ls_rip(:,ls_idx,ls_idx));
%         temporalBias_post{ii} = nan(size(ccg_ls_rip,2),size(ccg_ls_rip,3));
        warpedCorr_post{ii} = nan(size(ccg_ls_rip,2),size(ccg_ls_rip,3),250);
        warpedCorr_post_shuf{ii} = nan(size(ccg_ls_rip,2),size(ccg_ls_rip,3),1,250);
        for i = 1:size(ccg_ls_rip,2)
            for j = 1:size(ccg_ls_rip,3)
                if sum(ccg_ls_rip(:,i,j)) > 5 & sh(i,j) == 0 & ls_idx(i) == 1 & ls_idx(j) == 1 & sum(ccg_beh(:,i,j)) > 5
               
                b(i,j) = (sum(ccg_beh(51:250,i,j)))-(sum(ccg_beh(252:451,i,j)));
                br(i,j) = (sum(ccg_ls_rip(51:250,i,j)))-(sum(ccg_ls_rip(252:451,i,j)));
                
                
               %% time warping replay and temporal bias
                temporalBias_post{ii}(i,j) = single(sum(ccg_ls_rip(51:250,i,j) - ccg_ls_rip(252:451,i,j))) ./single(sum(ccg_ls_rip(51:250,i,j) + ccg_ls_rip(252:451,i,j))); 
%                 if sum(ccg_beh(:,i,j)) > 5
                    c1 = ccg_beh(:,i,j);
                    c2 = ccg_ls_rip(:,i,j);
                    for w=132:250
    %                    compression(w) = 240/(length(w:501-w+1));
                       warpedCorr_post{ii}(i,j,w) = single(corr((c1(131:371)),makeLength((c2(w:end-w+1)),241)'));
                       for iter = 1

                       warpedCorr_post_shuf{ii}(i,j,iter,w) = single(corr((c1(131:371)),circshift(makeLength((c2(w:end-w+1)),241),randi(240))'));

                       end
                    end
%                 end

                sma = Smooth(ccg_ls_rip(:,i,j),5);
                sma = sma ./ sum(sma);
                la = Smooth(ccg_ls_rip(:,i,j),200);
                la = la ./ sum(la);
                excess_post{ii}(i,j) = sma(251)-la(251);
                excess_max_post{ii}(i,j) = max(sma(225:275))-max(la(225:275));
                
                % jitter vers
                sma = Smooth(ccg_ls_rip_jitter(:,i,j),5);
                sma = sma ./ sum(sma);
                la = Smooth(ccg_ls_rip_jitter(:,i,j),200);
                la = la ./ sum(la);
                excess_post_jitter{ii}(i,j) = sma(251)-la(251);
                
                
                        
                nSpikes_post{ii}(i,j) = sum(ccg_ls_rip(:,i,j));
                else
                    warpedCorr_post{ii}(i,j,1:250) = single(NaN);
                    warpedCorr_post_shuf{ii}(i,j,1,1:250) = single(NaN);
                    temporalBias_post{ii}(i,j) = NaN;
                    b(i,j) = NaN;
                    br(i,j) = NaN;
                    excess_post{ii}(i,j) = NaN;
                    excess_post_jitter{ii}(i,j) = NaN;
                    excess_max_post{ii}(i,j) = NaN;
                    nSpikes_post{ii}(i,j) = 0;
                end
                
                
            end
        end
%         subplot(2,2,1)
    %     scatter(b(~sh(1:50,1:50)),br(~sh(1:50,1:50)))
        post_bias(ii) = (corr(b(~sh(ls_idx,ls_idx)),br(~sh(ls_idx,ls_idx)),'rows','complete'));
        post_bias_N(ii) = sum(~isnan(br(~sh(ls_idx,ls_idx))) & ~isnan(b(~sh(ls_idx,ls_idx))));
    end
    end
    clear b br  
    
    
    nEvtsPre(ii) = length(idx);

    idx = find(ripples.peaks<behavior.events.trialIntervals(1,1));
    nEvtsPost(ii) = length(idx);
   if ~isempty(idx)
    for i=1:length(spikes.times)
        spk.times{i} = Restrict(spikes.times{i},[ripples.timestamps(idx,1)-.05 ripples.timestamps(idx,2)+.05]);
%         spk.times{i} = [Restrict(spikes.times{i},[0 ripples.timestamps(idx(end),2)+.05])];
%         spk.times{i} = [Restrict(spikes.times{i},[behavior.events.trialIntervals(1,1) behavior.events.trialIntervals(end)])];
    end
    [times groups] = spikes2sorted(spk.times); clear spk
    if ~isempty(times) & sum(sum(~sh(ls_idx,ls_idx))) > 0
        [ccg_ls_rip t] = CCG(times,groups,'binSize',.001,'duration',.5,'norm','rate');
        ccg_ls_rip(end,length(spikes.times),length(spikes.times))=0;
        
        times_jitter = times + rand(length(times),1)./5;
        [ccg_ls_rip_jitter t] = CCG(times_jitter,groups,'binSize',.001,'duration',.5,'norm','rate');
        ccg_ls_rip_jitter(end,length(spikes.times),length(spikes.times))=0;
        
        ccg_pre{ii} = single(ccg_ls_rip(:,ls_idx,ls_idx));
%         temporalBias_pre{ii} = nan(size(ccg_ls_rip,2),size(ccg_ls_rip,3));
        warpedCorr_pre{ii} = nan(size(ccg_ls_rip,2),size(ccg_ls_rip,3),250);
        warpedCorr_pre_shuf{ii} = nan(size(ccg_ls_rip,2),size(ccg_ls_rip,3),1,250);
        for i = 1:size(ccg_ls_rip,2)
            for j = 1:size(ccg_ls_rip,3)
                if sum(ccg_ls_rip(:,i,j)) > 5 & sh(i,j) == 0 & ls_idx(i) == 1 & ls_idx(j) == 1 & sum(ccg_beh(:,i,j)) > 5
                    
                b(i,j) = (sum(ccg_beh(51:250,i,j)))-(sum(ccg_beh(252:451,i,j)));
                br(i,j) = (sum(ccg_ls_rip(51:250,i,j)))-(sum(ccg_ls_rip(252:451,i,j)));

                 %% time warping replay and temporal bias
                    
                temporalBias_pre{ii}(i,j) = single(sum(ccg_ls_rip(51:250,i,j) - ccg_ls_rip(252:451,i,j))) ./ single(sum(ccg_ls_rip(51:250,i,j) + ccg_ls_rip(252:451,i,j))); 
%                 if sum(ccg_beh(:,i,j)) > 5
                    c1 = ccg_beh(:,i,j);
                    c2 = ccg_ls_rip(:,i,j);
                    for w=132:250
                       warpedCorr_pre{ii}(i,j,w) = single(corr((c1(131:371)),makeLength((c2(w:end-w+1)),241)'));
                       for iter = 1
                           warpedCorr_pre_shuf{ii}(i,j,iter,w) = single(corr((c1(131:371)),circshift(makeLength((c2(w:end-w+1)),241),randi(240))'));
                       end
                    end
%                 end
            
                sma = Smooth(ccg_ls_rip(:,i,j),5);
                sma = sma ./ sum(sma);
                la = Smooth(ccg_ls_rip(:,i,j),200);
                la = la ./ sum(la);
                excess_pre{ii}(i,j) = sma(251)-la(251);
                excess_max_pre{ii}(i,j) = max(sma(225:275))-max(la(225:275));
                
                % jitter vers
                sma = Smooth(ccg_ls_rip_jitter(:,i,j),5);
                sma = sma ./ sum(sma);
                la = Smooth(ccg_ls_rip_jitter(:,i,j),200);
                la = la ./ sum(la);
                excess_pre_jitter{ii}(i,j) = sma(251)-la(251);
                
                
                
                        
                nSpikes_pre{ii}(i,j) = sum(ccg_ls_rip(:,i,j));
                else
                    nSpikes_pre{ii}(i,j) = 0;
                    warpedCorr_pre{ii}(i,j,1:250) = single(NaN);
                    warpedCorr_pre_shuf{ii}(i,j,1,1:250) = single(NaN);
                    temporalBias_pre{ii}(i,j) = NaN;
                    b(i,j) = NaN;
                    br(i,j) = NaN;
                    excess_pre{ii}(i,j) = NaN;
                    excess_pre_jitter{ii}(i,j) = NaN;
                    excess_max_pre{ii}(i,j) = NaN;
                end
                
               
                
            end
        end
%         subplot(2,2,2)
    %     scatter(b(~sh(1:50,1:50)),br(~sh(1:50,1:50)))
        pre_bias(ii) = (corr(b(~sh(ls_idx,ls_idx)),br(~sh(ls_idx,ls_idx)),'rows','complete'));
        pre_bias_N(ii) = sum(~isnan(br(~sh(ls_idx,ls_idx))) & ~isnan(b(~sh(ls_idx,ls_idx))));
    end
   end
    
    

    if ~isempty(excess_post{ii}) & ~isempty(excess_pre{ii})
        pre_excess_corr = [pre_excess_corr; (excess_pre{ii}(~sh))];
        post_excess_corr = [post_excess_corr; (excess_post{ii}(~sh))];
        beh_excess_corr = [beh_excess_corr; (excess_beh{ii}(~sh))];
        
        pre_excess_jitter_corr = [pre_excess_jitter_corr; (excess_pre_jitter{ii}(~sh))];
        post_excess_jitter_corr = [post_excess_jitter_corr; (excess_post_jitter{ii}(~sh))];
        beh_excess_jitter_corr = [beh_excess_jitter_corr; (excess_beh_jitter{ii}(~sh))];
        
        fieldCounts = [fieldCounts; fieldComp{ii}(~sh)];
    end
    end
    clear b br sh spk
    end
    ii
%     cd D:\datasets\ripples_LS
    cd /home/david/datasets/ripples_LS
% save('C:\Users\SB13FLLT001\Dropbox\data\HPC_replay_dataset_unsmoothed.mat','-v7.3')
end


idx = find(post_excess_corr==0 ...
           | pre_excess_corr==0 | ...
           isnan(post_excess_corr) ...
           | isnan(pre_excess_corr) ...
           | isnan(beh_excess_corr) ...
           | beh_excess_corr==0);
post_excess_corr(idx)=nan;
pre_excess_corr(idx)=nan;
beh_excess_corr(idx)=nan;


save('/home/david/Dropbox\data\HPC_replay_dataset_unsmoothed_5spks_allThresh_20190826.mat','-v7.3')

