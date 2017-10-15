function [] = scratch_10_11_17_compile_ls_hpc_noise_corrs(RECORDING)

cd(RECORDING)
d = dir('*firingMaps.cellinfo.mat');
sessionInfo = bz_getSessionInfo;

load([d.name])
for c1 = 1:length(firingMaps.UID)
  for c2 = 1:length(firingMaps.UID)
      if strcmp(firingMaps.region{c1},'ls') & strcmp(firingMaps.region{c2},'hpc')
          for cond = 1:length(firingMaps.rateMaps)
              if sum(sum(firingMaps.countMaps{cond}(c1,:,:))) > 20 & sum(sum(firingMaps.countMaps{cond}(c2,:,:))) > 20 & size(firingMaps.rateMaps{cond},2) > 10
                  m1 = squeeze(mean(firingMaps.rateMaps{cond}(c1,:,:)));
                  m2 = squeeze(mean(firingMaps.rateMaps{cond}(c2,:,:)));
                  for trial = 1:size(firingMaps.rateMaps{cond},2)
                     t1(trial,:) = [squeeze(firingMaps.rateMaps{cond}(c1,trial,:))-m1]; 
                     t2(trial,:) = [squeeze(firingMaps.rateMaps{cond}(c2,trial,:))-m2]; 
                  end
                  for bin = 1:size(t1,2)
                     noiseCorr(c1,c2,cond,bin) = corr(t1(:,bin),t2(:,bin),'rows','complete'); 
                     for iter = 1:100
                         r = randperm(trial);
                         noiseCorr_shuffle(c1,c2,cond,iter,bin) = corr(t1(r,bin),t2(:,bin),'rows','complete'); 
                     end
                  end
                  clear t1 t2 m1 m2
              end
          end             
      end         
  end
end
    
save([RECORDING '/' sessionInfo.FileName '.noiseCorrs.mat'],'noiseCorr','noiseCorr_shuffle')
   