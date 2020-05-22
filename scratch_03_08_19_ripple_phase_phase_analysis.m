addpath(genpath('/home/david/Dropbox/code'))
cd /home/david/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
clear all
ripReward = nan(98,201);

d = dir('*201*');

for rec = 1:length(d)
    cd(d(rec).name)
    coherence_trig{rec} = [];
    if exist(([d(rec).name '.behavior.mat']))
        sessionInfo = bz_getSessionInfo;
        sessionInfo.FileName;
        animal{rec} = sessionInfo.animal;
        if ~isempty(sessionInfo.ca1) && ~isempty(sessionInfo.ls)
        ca1_ripples = bz_LoadEvents(pwd,'CA1Ripples');
        ls_ripples = bz_LoadEvents(pwd,'LSRipples');
        ca1 = bz_GetLFP(sessionInfo.ca1);
        latS = bz_GetLFP(sessionInfo.ls);
        if ~isempty(ca1) && ~isempty(latS) && ~isempty(ca1_ripples) && ~isempty(ls_ripples)
            
            [coh phase t f] =  MTCoherogram([ca1.timestamps sin(angle(hilbert(bz_Filter(double(ca1.data),'passband',[120 200],'filter','butter'))))],...
                                            [latS.timestamps sin(angle(hilbert(bz_Filter(double(latS.data),'passband',[120 200],'filter','butter'))))]...
                                            ,'range',[120 200],'window',.05);
ct = []; 
ctc = []; 
cts = [];
ctcs = [];
            
        parfor rip = 1:length(ca1_ripples.timestamps)
%            idx = find(InIntervals(ls_ripples.peaks,[ca1_ripples.peaks(rip)-.025 ca1_ripples.peaks(rip)+.025]));
%             if ~isempty(idx)

            
                [a b] = min(abs(ca1_ripples.peaks(rip)-ca1.timestamps));
                start = b - 25;
                stop = b + 25;
                [a loc] = min(abs(t-ca1_ripples.timestamps(rip)));
                
                if start > 0 & stop < length(ca1.data) & loc > 500 & loc < length(coh)-500
%                     phase_ca1 = angle(hilbert(bz_Filter(double(ca1.data(start:stop)),'passband',[120 200],'filter','butter')));
%                     phase_ls = angle(hilbert(bz_Filter(double(latS.data(start:stop)),'passband',[120 200],'filter','butter')));
%                     ccg{rec}(rip,:) = crosscorr(sin(phase_ca1),sin(phase_ls),25);

                    ct(rip,:) = mean(coh(:,loc-100:loc+100),1);
                    loc_shuf = loc + randi(100)-50;
                    cts(rip,:) = mean(coh(:,loc_shuf-100:loc_shuf+100),1);
%                     for iter = 1:100
%                          ccg_shuf{rec}(rip,iter,:) = crosscorr(circshift(sin(phase_ca1),randi(51)),sin(phase_ls),25);
%                     end
                end
%             end

%% only coupled events now
            idx = find(InIntervals(ls_ripples.peaks,[ca1_ripples.peaks(rip)-.025 ca1_ripples.peaks(rip)+.025]));
            if ~isempty(idx)
                [a b] = min(abs(ca1_ripples.peaks(rip)-ca1.timestamps));
                start = b - 25;
                stop = b + 25;
                [a loc] = min(abs(t-ca1_ripples.timestamps(rip)));
                
                if start > 0 & stop < length(ca1.data) & loc > 100 & loc < length(coh)-100
%                     phase_ca1 = angle(hilbert(bz_Filter(double(ca1.data(start:stop)),'passband',[120 200],'filter','butter')));
%                     phase_ls = angle(hilbert(bz_Filter(double(latS.data(start:stop)),'passband',[120 200],'filter','butter')));
%                     ccg_coupled{rec}(rip,:) = crosscorr(sin(phase_ca1),sin(phase_ls),25);

                    ctc(rip,:) = mean(coh(:,loc-100:loc+100),1);
                    loc_shuf = loc + randi(100)-50;
                    ctcs(rip,:) = mean(coh(:,loc_shuf-100:loc_shuf+100),1);
%                     for iter = 1:100
%                          ccg_shuf_coupled{rec}(rip,iter,:) = crosscorr(circshift(sin(phase_ca1),randi(51)),sin(phase_ls),25);
%                     end
                end
            end
        end
        coherence_trig{rec} = ct;
        coherence_trig_shuf{rec} = cts;
        coherence_trig_coupled{rec} = ctc;
        coherence_trig_shuf_coupled{rec} = ctcs; clear ct*
        end
        end
        
    end
    plot(mean(cell2mat(coherence_trig')))
    pause(1)
    save('/home/david/Dropbox/data/phase_phase_coherence_10ms.mat','-v7.3')
   cd /home/david/datasets/ripples_LS 

end

