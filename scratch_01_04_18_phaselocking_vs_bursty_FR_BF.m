ranges = [4 9; 20 35; 40 60; 80 150];
load('rates_lfp_phase_lock_compiled_data_4_with_list_temp.mat', 'all_rates')
load('rates_lfp_phase_lock_compiled_data_4_with_list_temp.mat', 'raw_LFP')
load('rates_lfp_phase_lock_compiled_data_4_with_list_temp.mat', 'raw_spikes')
% load('bf_dataset_burstiness_FR_phaselocking.mat')
for cell = 1:780
    ntrials = size(raw_LFP{cell},1);
    lfp.data = reshape(raw_LFP{cell}(1:ntrials,:)',ntrials*6001,1);
    lfp.timestamps = [0:1/1000:length(lfp.data)/1000]';
    lfp_flipped = lfp;
    lfp_flipped.data = flipud(lfp_flipped.data);
    spks{1} = reshape(raw_spikes{cell}(1:ntrials,:)',ntrials*6001,1);
    spks{1} = find(spks{1})./1000;
    if spks{1}(1)==0
        spks{1} = spks{1}(2:end);
    end
   for freq = 1:4  
       resultants(cell,freq,1)=nan;
       if length(spks{1})>30
           iter = 1;
%        for iter = 1:10
%        r = randperm(length(spks{1}));
%        s{1} = spks{1}(r(1:300));
       s{1} = spks{1};
        phaselockingdata{cell}{freq} = bz_PhaseModulation(s,lfp,[ranges(freq,:)],...
            'samplingRate',1000,'powerThresh',0,'plotting',false,'saveMat',false);
         phaselockingdata_flipped{cell}{freq} = bz_PhaseModulation(s,lfp_flipped,[ranges(freq,:)],...
            'samplingRate',1000,'powerThresh',0,'plotting',false,'saveMat',false);
        resultants(cell,freq,iter) = phaselockingdata{cell}{freq}.phasestats.r;
        resultants_flipped(cell,freq,iter) = phaselockingdata_flipped{cell}{freq}.phasestats.r;
       end
%        end
   end
   FR(cell) = mean(spks{1}).*1000./length(spks{1});
   isi(cell,:) = hist(diff(spks{1}),0:.001:.25);
   burstiness(cell) = mean(isi(cell,1:25))./mean(isi(cell,25:150));
%    r = squeeze(mean(resultants,3));
%    rf = squeeze(mean(resultants_flipped,3));
%    
%    for freq = 1:4
%        subplot(2,2,freq)
%        scatter(mean(resultants(:,freq,:),3),log(FR),'.k')
%        [a b] = corr(r(:,i),FR','rows','complete');
%        title([a b])
%        [a b] = corr(rf(:,i),FR','rows','complete')
% %        [a b] = corr(r(:,i),burstiness','rows','complete')
%    end
%    pause(.01)
   cell
end
save('bf_dataset_burstiness_FR_phaselocking_non_normed.mat','-v7.3')



















