clf
hpc_phase = []; hpc_rate = [];
ls_rate = []; ls_phase = [];
d  = dir('*201*');
for i=1:length(d)
   cd(d(i).name) 
    if ~isempty(dir('*positionDecodingGLM*'))
        xml = LoadParameters;
        load([xml.FileName '.positionDecodingGLM.cellinfo.mat'])
        if isfield(positionDecodingGLM,'dateRun')
        conditions = length(unique(positionDecodingGLM.results{1}.condition));
        for cell =1:length(positionDecodingGLM.results)
            worstFit = max(positionDecodingGLM.results{cell}.mse_rate);
            for cond = 1:length(conditions)
               rows = find(positionDecodingGLM.results{cell}.condition==cond);
               [a b] =min(positionDecodingGLM.results{cell}.mse_phase_all(rows));
%                [a_sin b_sin] =min(positionDecodingGLM.results{cell}.mse_phase_sin(rows));
%                [a b] =min(positionDecodingGLM.results{cell}.mse_phase(rows));
%                bb = [b_cos b_sin b];
%                [a c] = min([a_cos a_sin a]);
%                b = bb(c);clear bb
               [aa bb] =min(positionDecodingGLM.results{cell}.mse_rate(rows));
               first500ms = find(ismember(positionDecodingGLM.results{cell}.tau(rows),1:500));
               if strcmp(positionDecodingGLM.region{cell},'hpc')
                   if positionDecodingGLM.results{cell}.mse_phase_all_pval(rows(b)) <.05 || ...
                           positionDecodingGLM.results{cell}.mse_rate(rows(bb)) <.05
                   subplot(2,4,1)
                   scatter(a./worstFit,aa./worstFit,'.k')
                   ylabel('rate')
                   xlabel('phase')
                   hold on
                   title('best fit (normed mse), any timescale')
                   axis([0 1 0 1])
                   
                   subplot(2,4,2)
                   scatter(positionDecodingGLM.results{cell}.tau(rows(b)),positionDecodingGLM.results{cell}.tau(rows(bb)),'.k')
                   hold on
                   title('best window (ms)')
                   ylabel('optimal rate time scale')
                   xlabel('optimal phase time scale')
                   
                   subplot(2,4,3)
                   min_mse_rate = min(positionDecodingGLM.results{cell}.mse_rate(rows(first500ms)));
                   min_mse_phase_all = min(positionDecodingGLM.results{cell}.mse_phase_all(rows(first500ms)));
                   hpc_phase=[hpc_phase;a./worstFit];
                   hpc_rate=[hpc_rate;aa./worstFit];
                   histogram(hpc_phase,100,'Normalization','pdf')
                   hold on
                   histogram(hpc_rate,100,'Normalization','pdf')
                   hold off
%                    scatter(min_mse_phase_all./worstFit,min_mse_rate./worstFit,'.k')
%                    hold on
%                    title('best fit (mse), 0-500 ms')
%                    xlabel('1-500 ms, phase normed MSE')
%                    ylabel('1-500 ms, rate normed MSE')
%                     axis([0 1 0 1])
                    
                   subplot(2,4,4)
                   plot(positionDecodingGLM.results{cell}.tau(rows(b)),a./worstFit,'.g')
                   hold on
                   plot(positionDecodingGLM.results{cell}.tau(rows(bb)),aa./worstFit,'.r')
                   title('tau vs mse trough')
                   end
               elseif strcmp(positionDecodingGLM.region{cell},'ls')    
                     if positionDecodingGLM.results{cell}.mse_phase_all_pval(rows(b)) <.05 || ...
                           positionDecodingGLM.results{cell}.mse_rate(rows(bb)) <.05
                   subplot(2,4,5)
                   scatter(a./worstFit,aa./worstFit,'.k')
                   ylabel('rate')
                   xlabel('phase')
                   hold on
                   title('best fit (normed mse), any timescale')
                   axis([0 1 0 1])
                   
                   subplot(2,4,6)
                   scatter(positionDecodingGLM.results{cell}.tau(rows(b)),positionDecodingGLM.results{cell}.tau(rows(bb)),'.k')
                   hold on
                   title('best window (ms)')
                   ylabel('optimal rate time scale')
                   xlabel('optimal phase time scale')
                   
                   subplot(2,4,7)
                   min_mse_rate = min(positionDecodingGLM.results{cell}.mse_rate(rows(first500ms)));
                   min_mse_phase_all = min(positionDecodingGLM.results{cell}.mse_phase_all(rows(first500ms)));
                   ls_phase=[ls_phase;a./worstFit];
                   ls_rate=[ls_rate;aa./worstFit];
                   histogram(ls_phase,100,'Normalization','pdf')
                   hold on
                   histogram(ls_rate,100,'Normalization','pdf')
                   hold off
%                    scatter(min_mse_phase_all./worstFit,min_mse_rate./worstFit,'.k')
%                    hold on
%                    title('best fit (mse), 0-500 ms')
%                    xlabel('1-500 ms, phase normed MSE')
%                    ylabel('1-500 ms, rate normed MSE')
%                     axis([0 1 0 1])
                    
                   subplot(2,4,8)
                   plot(positionDecodingGLM.results{cell}.tau(rows(b)),a./worstFit,'.g')
                   hold on
                   plot(positionDecodingGLM.results{cell}.tau(rows(bb)),aa./worstFit,'.r')
                   title('tau vs mse trough')
                   end
               end
            end       
        end
        end
    end
    pause(.1)
    cd /home/david/datasets/lsDataset/
end