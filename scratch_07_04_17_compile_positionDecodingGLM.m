clf
hpc_phase = []; hpc_rate = [];
ls_rate = []; ls_phase = [];
d  = dir('*201*');
for i=1:length(d)
   cd(d(i).name) 
    if ~isempty(dir('*positionDecodingGLM_binnedspace.cell*'))
        xml = LoadParameters;
        load([xml.FileName '.positionDecodingGLM_binnedspace.cellinfo.mat'])
        positionDecodingGLM=positionDecodingGLM_binnedspace;
        if isfield(positionDecodingGLM,'dateRun')
        conditions = length(unique(positionDecodingGLM.results{1}.condition));
        for cell =1:length(positionDecodingGLM.results)
            t_rate = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_rate',...
                'GroupingVariables',{'tau','condition'});
            t_phase = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_phase_all',...
                'GroupingVariables',{'tau','condition'});
            tab = join(t_rate,t_phase);
            for cond = 1:conditions
               rows = find(tab.condition==cond);
               
%                worstFit = max([positionDecodingGLM.results{cell}.mse_rate(rows)]);
               
               [a b] =min(t_phase.mean_mse_phase_all(rows));
%                [a_sin b_sin] =min(positionDecodingGLM.results{cell}.mse_phase_sin(rows));
%                [a b] =min(positionDecodingGLM.results{cell}.mse_phase(rows));
%                bb = [b_cos b_sin b];
%                [a c] = min([a_cos a_sin a]);
%                b = bb(c);clear bb
               [aa bb] =min(tab.mean_mse_rate(rows));
               
               first500ms = find(ismember(tab.tau(rows),1:500));
               
               min_mse_rate = min(tab.mean_mse_phase_all(rows(first500ms)));
               min_mse_phase_all = min(tab.mean_mse_rate(rows(first500ms)));
               
%                if b ~= length(rows) & bb ~= length(rows) & b ~= 1 & bb ~= 1
               if strcmp(positionDecodingGLM.region{cell},'hpc')
%                    if positionDecodingGLM.results{cell}.mse_phase_all_pval(rows(b)) <.05 || ...
%                            positionDecodingGLM.results{cell}.mse_rate_pval(rows(bb)) <.05
                   
                   subplot(2,4,1)
                   scatter(min_mse_phase_all,min_mse_rate,'.k')
                   ylabel('rate')
                   xlabel('phase')
                   hold on
                   title('best fit (normed mse), any timescale')
%                    axis([0 1 0 1])
                   
                   subplot(2,4,2)
                   scatter(tab.tau(rows(b)),tab.tau(rows(bb)),'.k')
                   hold on
                   title('best window (ms)')
                   ylabel('optimal rate time scale')
                   xlabel('optimal phase time scale')
                   
                   subplot(2,4,3)
                   hpc_phase=[hpc_phase;min_mse_phase_all];
                   hpc_rate=[hpc_rate;min_mse_rate];
                   histogram(hpc_phase,100,'Normalization','pdf','FaceColor','g'); .4:.05:1;
                   hold on
                   histogram(hpc_rate,100,'Normalization','pdf','FaceColor','r')
                   set(gca,'yscale','log')
                   hold off
%                    scatter(min_mse_phase_all./worstFit,min_mse_rate./worstFit,'.k')
%                    hold on
%                    title('best fit (mse), 0-500 ms')
%                    xlabel('1-500 ms, phase normed MSE')
%                    ylabel('1-500 ms, rate normed MSE')
%                     axis([0 1 0 1])
                    
                   subplot(2,4,4)
                   plot(tab.tau(rows(b)),min_mse_phase_all,'.g')
                   hold on
                   plot(tab.tau(rows(bb)),min_mse_rate,'.r')
                   title('tau vs mse trough')
%                    end
               elseif strcmp(positionDecodingGLM.region{cell},'ls')    
%                      if positionDecodingGLM.results{cell}.mse_phase_all_pval(rows(b)) <.05 || ...
%                            positionDecodingGLM.results{cell}.mse_rate(rows(bb)) <.05
                   subplot(2,4,5)
                   scatter(min_mse_phase_all,min_mse_rate,'.k')
                   ylabel('rate')
                   xlabel('phase')
                   hold on
                   title('best fit (normed mse), any timescale')
%                    axis([0 1 0 1])
                   

                   subplot(2,4,7)
                   ls_phase=[ls_phase;min_mse_phase_all];
                   ls_rate=[ls_rate;min_mse_rate];
                   histogram(ls_phase,100,'Normalization','pdf','FaceColor','g')
                   hold on
                   histogram(ls_rate,100,'Normalization','pdf','FaceColor','r')
                   set(gca,'yscale','log')
                   hold off
%                    scatter(min_mse_phase_all./worstFit,min_mse_rate./worstFit,'.k')
%                    hold on
%                    title('best fit (mse), 0-500 ms')
%                    xlabel('1-500 ms, phase normed MSE')
%                    ylabel('1-500 ms, rate normed MSE')
%                     axis([0 1 0 1])
                    
                   subplot(2,4,6)
                   histogram((hpc_phase-hpc_rate)./(hpc_phase+hpc_rate),[-.4:.01:.4],'Normalization','pdf','FaceColor','r')
                   set(gca,'yscale','log')
                   hold on
                   histogram((ls_phase-ls_rate)./(ls_phase+ls_rate),[-.4:.01:.4],'Normalization','pdf','FaceColor','g')
                   hold off
                   subplot(2,4,8)
                   plot(tab.tau(rows(b)),min_mse_phase_all,'.g')
                   hold on
                   plot(tab.tau(rows(bb)),min_mse_rate,'.r')
                   title('tau vs mse trough')
%                      end
               end
%                end   
            pause(.1)
            end
        end
        end
    end
    
    cd /home/david/datasets/lsDataset/
end