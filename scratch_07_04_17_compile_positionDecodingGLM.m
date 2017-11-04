clf
hpc_phase = []; hpc_rate = [];
ls_rate = []; ls_phase = [];
ls_tau_phase = []; ls_tau_rate = [];
hpc_tau_phase =  []; hpc_tau_rate = [];
hpc_phase_pval = [];
hpc_rate_pval = [];
ls_phase_pval = [];
ls_rate_pval = [];
chance =[];
ls_rec =[]; hpc_rec = [];
hpc_an = []; ls_an=[];
ls_depth=[];
hpc_depth=[];
ls_phase_info = []; hpc_phase_info =[];
ls_rate_info =[]; hpc_rate_info = [];
hpc_cell=[];
ls_cell=[];
hpc_shank = []; ls_shank = [];
hpc_field = [];
                   
hpc_mean_phase = zeros(101,1);
ls_mean_phase = zeros(101,1);
hpc_mean_rate = zeros(101,1);
ls_mean_rate = zeros(101,1);
d  = dir('*201*');
for i=1:length(d)
    i
   cd(d(i).name) 
   sessionInfo = bz_getSessionInfo;
   animal = sessionInfo.animal;
%    animal = strsplit(animal,'/');
%    animal = animal{end-1};
% animal = 1;
    if ~isempty(dir('*positionDecodingGLM_binnedspace_box.cell*')) & exist([d(i).name '.placeFields.10_pctThresh.mat'])
        sessionInfo = bz_getSessionInfo;
%         load([d(i).name '.firingMaps.cellinfo.mat'],'firingMaps') 
        load([d(i).name '.placeFields.10_pctThresh.mat'],'fields') 
        spikes = bz_GetSpikes;
%         load([d(i).name '.olypherInfo.cellinfo.mat'],'olypherInfo') 
        b = dir('*.behavior.mat');
        load(b(1).name);
        nBins = round(length(behavior.events.map{1}.x));
        load([sessionInfo.FileName '.positionDecodingGLM_binnedspace_box.cellinfo.mat'])
        positionDecodingGLM=positionDecodingGLM_binnedspace_box;
        if isfield(positionDecodingGLM,'dateRun')
        conditions = length(unique(behavior.events.trialConditions));
        for cell =1:length(positionDecodingGLM.results)
            if ~isempty(positionDecodingGLM.results{cell})
        t_rate = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_rate',...
            'GroupingVariables',{'tau','condition'});
        t_phase = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_phase_all',...
            'GroupingVariables',{'tau','condition'});
        t_chance = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_chance',...
            'GroupingVariables',{'tau','condition'});
        tab = join(join(t_rate,t_phase),t_chance);
        
        t_phase_pval = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_phase_all_pval',...
            'GroupingVariables',{'tau','condition'});
        t_rate_pval = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_rate_pval',...
            'GroupingVariables',{'tau','condition'});       
        pvals = join(t_phase_pval,t_rate_pval);
        
            for cond = 1:conditions
                if sum(behavior.events.trialConditions==cond) >= 12 %%%%%%%%%%%%%%%%%%%%%%%%%%
                    nTrials = sum(behavior.events.trialConditions==cond);
                    %% information theory stuff here
%                     rows = find(olypherInfo.results{cell}.condition==cond);
%                     cols = find(olypherInfo.results{cell}.discBins==4);
%                     rows = intersect(rows,cols);
%                     maxphaseInfo = max(olypherInfo.results{cell}.phasePeakInfo(rows)./nTrials);
%                     maxrateInfo = max(olypherInfo.results{cell}.ratePeakInfo(rows)./nTrials);
                    
                    
                    
                    %% carry on..
               rows = find(tab.condition==cond);
%                if sqrt(tab.mean_mse_phase_all(rows(1)))./nBins < .3 
                [a b] =min(tab.mean_mse_phase_all(rows));
                [aa bb] =min(tab.mean_mse_rate(rows));
                
%                 rows = intersect(rows,find(tab.tau==60));
                
                first500ms = find(ismember(tab.tau(rows),4:7));

%                 min_mse_rate = (min(tab.mean_mse_rate(rows(first500ms)))./mean(tab.mean_mse_chance(rows(first500ms))));
%                 min_mse_phase_all = (min(tab.mean_mse_phase_all(rows(first500ms)))./mean(tab.mean_mse_chance(rows(first500ms))));
                
                [min_mse_rate] = (mean(tab.mean_mse_rate(rows(first500ms)))-mean(tab.mean_mse_chance(rows)));
                [min_mse_phase_all] = (mean(tab.mean_mse_phase_all(rows(first500ms)))-mean(tab.mean_mse_chance(rows)));

                min_mse_chance = (mean(tab.mean_mse_chance(first500ms)))-mean(tab.mean_mse_chance(rows));

%                 max_mse_rate = sqrt(max(tab.mean_mse_rate(rows(first500ms))))./nBins;
%                 max_mse_phase_all = sqrt(max(tab.mean_mse_phase_all(rows(first500ms))))./nBins;
%                 if isempty(max_mse_phase_all)
%                     error();
%                 end
                if mean(tab.mean_mse_chance(rows)) > 3800
                   min_mse_rate = nan;
                   min_mse_phase_all = nan;
                   min_mse_chance = nan;
                end
%                if min_mse_phase_all < .33 & min_mse_rate < .33
%                if b ~= length(rows) & bb ~= length(rows) & b ~= 1 & bb ~= 1
                chance = [chance; min_mse_chance];
               if strcmp(positionDecodingGLM.region{cell},'hpc') | strcmp(positionDecodingGLM.region{cell},'ca3')  | strcmp(positionDecodingGLM.region{cell},'ca1') 
%                    if positionDecodingGLM.results{cell}.mse_phase_all_pval(rows(b)) <.05 || ...
%                            positionDecodingGLM.results{cell}.mse_rate_pval(rows(bb)) <.05
                   
%                    subplot(2,5,1)
%                    scatter(min_mse_phase_all,min_mse_rate,'.k')
%                    ylabel('rate')
%                    xlabel('phase')
%                    hold on
%                    title('best fit (normed mse), any timescale')
%                    axis([0 1 0 1])
%                    
%                    subplot(2,5,2)
%                    scatter(tab.tau(rows(b))./nBins,tab.tau(rows(bb))./nBins,'.k')
%                    hold on
%                    title('best window (ms)')
%                    ylabel('optimal rate time scale')
%                    xlabel('optimal phase time scale')
%                    
%                    subplot(2,5,3)
% hpc_phase_info = [hpc_phase_info; maxphaseInfo];
% hpc_rate_info = [hpc_rate_info; maxrateInfo];
                   hpc_phase=[hpc_phase;min_mse_phase_all];
                   hpc_rate=[hpc_rate;min_mse_rate];
                   hpc_tau_phase = [hpc_tau_phase;tab.tau(rows(b))./nBins];
                   hpc_tau_rate = [hpc_tau_rate;tab.tau(rows(bb))./nBins];
                   hpc_phase_pval = [hpc_phase_pval;pvals.mean_mse_phase_all_pval(rows(b),:)];
                   hpc_rate_pval = [hpc_rate_pval;pvals.mean_mse_rate_pval(rows(bb),:)];
                   hpc_rec = [hpc_rec; i];
                   hpc_an = [hpc_an; sum(double(animal))];
                   additionalDepth = find(sessionInfo.spikeGroups.groups{spikes.shankID(cell)}==spikes.maxWaveformCh(cell))*10;
                   hpc_depth = [hpc_depth;str2num(sessionInfo.depth)+additionalDepth];
                   hpc_cell = [hpc_cell; cell];
                   hpc_shank = [hpc_shank;spikes.shankID(cell)];
                   
                   if ~isempty(fields{cond}{cell}) & length(fields{cond}{cell}) == 1
                       hpc_field = [hpc_field;fields{cond}{cell}{1}.COM];
                   else
                       hpc_field = [hpc_field;nan];
                   end
%                    histogram(hpc_phase,0:.01:1,'Normalization','pdf','FaceColor','g'); .4:.05:1;
%                    hold on
%                    histogram(hpc_rate,0:.01:1,'Normalization','pdf','FaceColor','r')
%                    set(gca,'yscale','log')
%                    hold off
%                    
%                    subplot(2,5,4)
                   
%                    scatter(sqrt(a)./nBins,sqrt(aa)./nBins,'.k')
%                    hold on
%                    title('best fit (mse), 0-500 ms')
%                    xlabel('1-500 ms, phase normed MSE')
%                    ylabel('1-500 ms, rate normed MSE')
                    
%                    subplot(2,5,5)
%                    plot(tab.tau(rows(b))./nBins,min_mse_phase_all,'.g')
%                    hold on
%                    plot(tab.tau(rows(bb))./nBins,min_mse_rate,'.r')
%                    axis([0 .5 0 1])
%                    hpc_mean_phase = hpc_mean_phase + sqrt(tab.mean_mse_phase_all(rows))./nBins;
%                    plot(tab.tau(rows),hpc_mean_phase,'g')
%                    hpc_mean_rate = hpc_mean_rate + sqrt(tab.mean_mse_rate(rows))./nBins;
%                    hold on
%                    plot(tab.tau(rows),hpc_mean_rate,'r')
%                    hold off
                   
%                    title('tau vs mse trough')
%                    end
               elseif strcmp(positionDecodingGLM.region{cell},'ls') 
%                      if positionDecodingGLM.results{cell}.mse_phase_all_pval(rows(b)) <.05 || ...
%                            positionDecodingGLM.results{cell}.mse_rate(rows(bb)) <.05
%                    subplot(2,5,6)
%                    scatter(min_mse_phase_all,min_mse_rate,'.k')
%                    ylabel('rate')
%                    xlabel('phase')
%                    hold on
%                    title('best fit (normed mse), any timescale')
%                    axis([0 1 0 1])
%                    
% 
%                    subplot(2,5,7)
%                    scatter(tab.tau(rows(b))./nBins,tab.tau(rows(bb))./nBins,'.k')
%                    hold on
%                    title('best window (ms)')
%                    ylabel('optimal rate time scale')
%                    xlabel('optimal phase time scale')
%                    
%                    subplot(2,5,8)
% ls_phase_info = [ls_phase_info; maxphaseInfo];
% ls_rate_info = [ls_rate_info; maxrateInfo];
                   ls_phase=[ls_phase;min_mse_phase_all];
                   ls_rate=[ls_rate;min_mse_rate];
                   ls_tau_phase = [ls_tau_phase;tab.tau(rows(b))./nBins];
                   ls_tau_rate = [ls_tau_rate;tab.tau(rows(bb))./nBins];
                   ls_phase_pval = [ls_phase_pval;pvals.mean_mse_phase_all_pval(rows(b),:)];
                   ls_rate_pval = [ls_rate_pval;pvals.mean_mse_rate_pval(rows(bb),:)];
                   ls_rec = [ls_rec; i];
                   ls_an = [ls_an; sum(double(animal))];
                   additionalDepth = find(sessionInfo.spikeGroups.groups{spikes.shankID(cell)}==spikes.maxWaveformCh(cell))*10;
                   ls_depth = [ls_depth; str2num(sessionInfo.depth)+additionalDepth];
                   ls_cell = [ls_cell; cell];
                   ls_shank = [ls_shank;spikes.shankID(cell)];
%                    histogram(ls_phase,0:.01:1,'Normalization','pdf','FaceColor','g')
%                    hold on
%                    histogram(ls_rate,0:.01:1,'Normalization','pdf','FaceColor','r')
%                    set(gca,'yscale','log')
%                    hold off
%                    
%                    subplot(2,5,9)
% %                    scatter(sqrt(a)./nBins,sqrt(aa)./nBins,'.k')
% %                    hold on
% %                    title('best fit (mse), 0-500 ms')
% %                    xlabel('1-500 ms, phase normed MSE')
% %                    ylabel('1-500 ms, rate normed MSE')
% 
%                     
% %                    subplot(2,5,9)
% %                    histogram((hpc_phase-hpc_rate)./(hpc_phase+hpc_rate),[-.4:.01:.4],'Normalization','pdf','FaceColor','r')
% %                    set(gca,'yscale','log')
% %                    hold on
% %                    histogram((ls_phase-ls_rate)./(ls_phase+ls_rate),[-.4:.01:.4],'Normalization','pdf','FaceColor','g')
% %                    hold off
%                    
%                    subplot(2,5,10)
%                    plot(tab.tau(rows(b))./nBins,min_mse_phase_all,'.g')
%                    hold on
%                    plot(tab.tau(rows(bb))./nBins,min_mse_rate,'.r')
%                    axis([0 .5 0 1])
% % ls_mean_phase = ls_mean_phase + sqrt(tab.mean_mse_phase_all(rows))./nBins;
% %                     plot(tab.tau(rows),ls_mean_phase,'g')
% %                     hold on
% %                     ls_mean_rate = ls_mean_rate + sqrt(tab.mean_mse_rate(rows))./nBins;
% %                     plot(tab.tau(rows),ls_mean_rate,'r')
% %                     hold off
%                     title('tau vs mse trough')
%                      end
               end
%                end   
%                end
%                 end
%             pause(.01)
                end
            end
            end
        end
        end
    end
%     pause
%     cd('D:\Dropbox\datasets\lsDataset\')
cd('/home/david/datasets/lsDataset')
end

for i=1:length(d)
f = find(ls_rec==i);
ff = find(ls_phase(f)<1);
if ~isempty(ff)
anim(i) = ls_an(f(ff(1)));
depths(i) = ls_depth(f(ff(1)));
end
end
u = unique(anim);
for i=1:length(d)
f = find(hpc_rec==i);
ff = find(hpc_phase(f)<1);
hp(i) = mean(hpc_phase(f(ff)));
hr(i) = mean(hpc_rate(f(ff)));
f = find(ls_rec==i);
ff = find(ls_phase(f)<1);
if ~isempty(f(ff))
lmp(i) = nanmedian(ls_phase(f(ff)));
lmr(i) = nanmedian(ls_rate(f(ff)));
lr(i) = nanmean(ls_rate(f(ff)));
lp(i) = nanmean(ls_phase(f(ff)));
else
lmp(i) = 0;
lmr(i) = 0;
lr(i) = 0;
lp(i) = 0; 
end
end
% offsets = [0 0 -3000 -800 0 0 0];
for i=1:length(u)
subplot(4,2,i);hold on
f = find(anim==u(i));
plot(depths(f),lmp(f)-lmr(f),'.k')
[a  b] = corr(depths(f)',lmp(f)'-lmr(f)');
[x y] = polyfit(depths(f),lmp(f)-lmr(f),1);
y1 = polyval(x,depths(f));
plot(depths(f),y1,'r')
title(a)
end


