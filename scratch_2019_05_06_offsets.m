cd ~/datasets/ripples_LS
% cd E:\datasets\ripples_LS
d = dir('*201*');
c=1;
cc = 1;
specs = nan(25000,625,100);
specs_hpc = zeros(625,100);
specs_ls = zeros(625,100);


for i=1:98
    cd(d(i).name)
    if exist([d(i).name '.CA1Ripples.events.mat']) & exist([d(i).name '.LSRipples.events.mat'])
        rips_ls = load([d(i).name '.LSRipples.events.mat']);
        rips_hpc = load([d(i).name '.CA1Ripples.events.mat']);
        zamps = zscore(rips_hpc.ripples.data.peakAmplitude);
        zfreqs = (rips_ls.ripples.data.peakFrequency);
        zdurs = (rips_ls.ripples.data.duration);
        
        sessionInfo = bz_getSessionInfo;

        
        for r = 1:length(rips_hpc.ripples.peaks)
           [a b] = min(abs(rips_ls.ripples.peaks-rips_hpc.ripples.peaks(r)));
           
           offsets(c) = rips_ls.ripples.peaks(b)-rips_hpc.ripples.peaks(r);
           durs(c) = zdurs(b); %rips_ls.ripples.data.duration(r);
           freqs(c) = zfreqs(b); %rips_ls.ripples.data.peakFrequency(r);
           
           amps_hpc(c) = zamps(r);
           
           durs_hpc(c) = rips_hpc.ripples.data.duration(r);
           freq_hpc(c) =rips_hpc.ripples.data.peakFrequency(r);
           
           if r > 1 & r < length(rips_hpc.ripples.peaks)
            iri_hpc(c,:) = (diff(rips_hpc.ripples.peaks(r-1:r+1)));
           else
            iri_hpc(c,:) = [nan nan];
           end
           
           rec(c) = i;
           
%             if ~isempty(sessionInfo.ca1)
%                lfp = bz_GetLFP(sessionInfo.ca1,'intervals',[rips_ls.ripples.peaks(r)-.25 rips_ls.ripples.peaks(r)+.25]);
%                spec = bz_WaveSpec(lfp,'frange',[1 300],'nfreqs',100,'ncyc',3,'space','lin');
% 
%                specs(c,:,:) = abs(spec.data);
%                specs_ls = specs_ls + abs(spec.data);;
%             else
%             specs(c,:,:) = nan;    
%             end
           c=1+c;

        end
%         for r = 1:length(rips_hpc.ripples.peaks)
%             [a b] = min(abs(rips_ls.ripples.peaks-rips_hpc.ripples.peaks(r)));
%             offsets_hpc(cc) = rips_ls.ripples.peaks(b)-rips_hpc.ripples.peaks(r);
%             durs_hpc(cc) = rips_hpc.ripples.data.duration(r);
%             freq_hpc(cc) =rips_hpc.ripples.data.peakFrequency(r);
%              
%             if ~isempty(sessionInfo.ca1)
% %                lfp = bz_GetLFP(sessionInfo.ca1,'intervals',[rips_hpc.ripples.peaks(r)-.25 rips_hpc.ripples.peaks(r)+.25]);
% %                spec = bz_WaveSpec(lfp,'frange',[1 300],'nfreqs',100,'ncyc',3,'space','lin');
% %                specs_hpc = specs_hpc + abs(spec.data);
%             else
%             
%             end
%            cc=1+cc;
%         end
%         rLS(round(rips_ls.ripples.peaks*1000)) = 1;
%         rHPC(round(rips_hpc.ripples.peaks*1000)) = 1;
%         rHPC(length(rLS)) = 0;
%         rLS(length(rHPC)) = 0;
%         rLS(length(rHPC)) = 0;
% 
% 
%         [f(i) cv(i)] = granger_cause(rHPC,rLS,.05,10);
%         N(i) = length(rips_ls.ripples.peaks);
%         
%         [pk(i) loc(i)] =max(crosscorr(rHPC,rLS,100));
%         
%         [times groups] = spikes2sorted({rips_hpc.ripples.peaks,rips_ls.ripples.peaks});
%         [ccg{i} t] = CCG(times,groups,'binSize',.001,'duration',.2);
        
        clear rip5*
    end
%     if i > 2
%     for k = 1:100
%         sl(:,k) = zscore(specs_ls(:,k));
%         sh(:,k) = zscore(specs_hpc(:,k));
%     end
%     subplot(2,2,1);
%     imagesc(spec.timestamps,spec.freqs,sl');
%     subplot(2,2,2)
%     imagesc(spec.timestamps,spec.freqs,sh');
%     pause(.1)
%     end
    i
    cd ~/datasets/ripples_LS
% cd E:\datasets\ripples_LS
end