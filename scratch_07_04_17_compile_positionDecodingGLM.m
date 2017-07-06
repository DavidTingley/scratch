clf

d  = dir('*201*');
for i=1:length(d)
   cd(d(i).name) 
    if ~isempty(dir('*positionDecodingGLM*'))
        xml = LoadParameters;
        load([xml.FileName '.positionDecodingGLM.cellinfo.mat'])
        conditions = length(unique(positionDecodingGLM.results{1}.condition));
        for cell =1:length(positionDecodingGLM.results)
            worstFit = max(positionDecodingGLM.results{cell}.mse_rate);
            for cond = 1:length(conditions)
               rows = find(positionDecodingGLM.results{cell}.condition==cond);
               [a_cos b_cos] =min(positionDecodingGLM.results{cell}.mse_phase_cos(rows));
               [a_sin b_sin] =min(positionDecodingGLM.results{cell}.mse_phase_sin(rows));
               [a b] =min(positionDecodingGLM.results{cell}.mse_phase(rows));
               bb = [b_cos b_sin b];
               [a c] = min([a_cos a_sin a]);
               b = bb(c);clear bb
               [aa bb] =min(positionDecodingGLM.results{cell}.mse_rate(rows));
               first500ms = find(ismember(positionDecodingGLM.results{cell}.tau(rows),1:500));
               if strcmp(positionDecodingGLM.region{cell},'hpc')
                   subplot(2,2,1)
                   scatter(a./worstFit,aa./worstFit,'.k')
                   ylabel('rate')
                   xlabel('phase')
                   hold on
                   title('best fit (normed mse), any timescale')
                   axis([0 1 0 1])
                   
                   subplot(2,2,2)
                   if positionDecodingGLM.results{cell}.mse_phase_cos_pval(rows(b)) <.05 || ...
                           positionDecodingGLM.results{cell}.mse_rate(rows(bb)) <.05
                   scatter(positionDecodingGLM.results{cell}.tau(rows(b)),positionDecodingGLM.results{cell}.tau(rows(bb)),'.k')
                  
                   end
                   hold on
                   title('best window (ms)')
                   ylabel('optimal rate time scale')
                   xlabel('optimal phase time scale')
                   
                   subplot(2,2,3)
                   min_mse_rate = min(positionDecodingGLM.results{cell}.mse_rate(rows(first500ms)));
                   min_mse_phase_cos = min(positionDecodingGLM.results{cell}.mse_phase_cos(rows(first500ms)));
                   scatter(min_mse_phase_cos./worstFit,min_mse_rate./worstFit,'.k')
                   hold on
                   title('best fit (mse), 0-500 ms')
                   xlabel('1-500 ms, phase normed MSE')
                   ylabel('1-500 ms, rate normed MSE')
                    axis([0 1 0 1])
                    
                   subplot(2,2,4)
                   plot(positionDecodingGLM.results{cell}.tau(rows(b)),a./worstFit,'.g')
                   hold on
                   plot(positionDecodingGLM.results{cell}.tau(rows(bb)),aa./worstFit,'.r')
                   title('tau vs mse trough')
                   
               elseif strcmp(positionDecodingGLM.region{cell},'ls')                   
               end
            end       
        end
    end
    pause(.1)
    cd /home/david/datasets/lsDataset/
end