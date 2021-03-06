cd /home/david/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
% clear all
d = dir('*201*');
count=1;


wakeCrossRegionCoupling = nan(length(d),1);
NREMCrossRegionCoupling = nan(length(d),1);
REMCrossRegionCoupling = nan(length(d),1);
wakeLSRippleRate = nan(length(d),1);
NREMLSRippleRate = nan(length(d),1);
REMLSRippleRate = nan(length(d),1);
wakeLSRippleCount = nan(length(d),1);
NREMLSRippleCount = nan(length(d),1);
REMLSRippleCount = nan(length(d),1);
wakeCA1RippleRate = nan(length(d),1);
NREMCA1RippleRate = nan(length(d),1);
REMCA1RippleRate = nan(length(d),1);

count = 1;
ratio = [];
reg =[];

for rec = 1:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    sessionInfo.FileName;
    ls_ripples = bz_LoadEvents(pwd,'LSRipples');
    ca1_ripples = bz_LoadEvents(pwd,'CA1Ripples');
    sleep = bz_LoadStates(pwd,'SleepState');
    spikes = bz_GetSpikes('noprompt',true);
%     if exist([sessionInfo.FileName '.placeFields.20_pctThresh.mat']) & exist([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
%     load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
%     load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
    if ~isempty(ls_ripples) 
        nEvents(rec) = length(ls_ripples.peaks);
        if ~ isempty(ca1_ripples)
%         if length(ca1_ripples.peaks) > 10
%     % get coupling scores
    for i=1:length(ls_ripples.timestamps) % ca1 as ref for now...
       if min(abs(ca1_ripples.peaks-ls_ripples.peaks(i))) < .025
            coupling(i) = 1;
       else
            coupling(i) = 0;
       end
    end
%        for spk = 1:length(spikes.times)
%           spkCount(spk) = length(Restrict(spikes.times{spk},[ca1_ripples.peaks-.2  ca1_ripples.peaks+.2]));
%           spkCount_shuffle(spk) = length(Restrict(spikes.times{spk},[ca1_ripples.peaks-.2  ca1_ripples.peaks+.2]+randi(30,length(ca1_ripples.peaks),1)));
%           region(spk) = sum(double(spikes.region{spk}));
%           
%           % get spatial info
%           maxTau = max(positionDecodingMaxCorr_binned_box_mean.results{spk}.tau);
%         rows = find(positionDecodingMaxCorr_binned_box_mean.results{spk}.tau==20);
%         r = positionDecodingMaxCorr_binned_box_mean.results{spk}.mse_rate(rows);
%         p(:,1) = positionDecodingMaxCorr_binned_box_mean.results{spk}.mse_phase_cos(rows);
%         p(:,2) = positionDecodingMaxCorr_binned_box_mean.results{spk}.mse_phase(rows);
%         p(:,3) = positionDecodingMaxCorr_binned_box_mean.results{spk}.mse_phase_all(rows);
%         [t loc]=min(mean(p));
%         p = p(:,loc);
%         c = positionDecodingMaxCorr_binned_box_mean.results{spk}.mse_chance_rate(rows);
%         phaseCode(spk) = mean(p);
%         phaseCodeShuf(spk) = mean(c);
%         rateCode(spk) = mean(r);
%        end
%        spkCounts{rec} = spkCount;
%        spkCounts_shuffle{rec} = spkCount_shuffle;
%        phCode{rec} = phaseCode; 
%        rCode{rec} = rateCode;
%        phCodeShuf{rec} = phaseCodeShuf;
%        regions{rec} = region; 
%        power{rec} = ca1_ripples.peakNormedPower;
%        coup{rec} = coupling; 
%     ratio = [ratio, ((spkCount)+1)./((spkCount_shuffle)+1)];
%     reg = [reg, region];
%     histogram(ratio(reg==315 | reg == 245),'normalization','pdf','facecolor','k')
%     hold on
%     histogram(ratio(reg==223),'normalization','pdf','facecolor','m')
%     histogram(ratio(reg==247),'normalization','pdf','facecolor','r')
%     set(gca,'yscale','log')
%     hold off
%     pause(.01)
%     clear spkCount spkCount_shuffle region phaseCode* rateCode c p r t loc coupling
%         end
%     end
    
%     for spk = 1:length(spikes.times)
%         for cond = 1:length(fields)
%            if ~isempty(fields{cond}{spk})
%                f(spk) = 1;
%            else
%                f(spk) = 0;
%            end
%         end
%         [times groups] = spikes2sorted({spikes.times{spk},ca1_ripples.peaks(coupling==0)});
%         [ccg t] = CCG(times,groups,'binSize',.001,'duration',.5);
%         ccg = ccg ./ sum(coupling==0);
%         [times groups] = spikes2sorted({spikes.times{spk},ca1_ripples.peaks(coupling==1)});
%         [ccg_coupled t] = CCG(times,groups,'binSize',.001,'duration',.5);
%         ccg_coupled = ccg_coupled ./ sum(coupling==1);
%         if f(spk) == 1 & ~strcmp(spikes.region{spk},'ls') & sum(coupling) > 15
%             plot(ccg(:,1,2),'k')
%             hold on
%             plot(ccg_coupled(:,1,2),'r')
%             pause(.001)
%             hold off
%             not_coup(count,:) = ccg(:,1,2);
%             coupl(count,:) = ccg_coupled(:,1,2);
%             count = 1+count;
%         end
%         clear f
        end
    end
    % get wake rate
    wakeTime = double(sum(diff(sleep.ints.WAKEstate')));
    if ~isempty(ls_ripples)
        wakeLSRippleRate(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.WAKEstate)))./wakeTime;
        wakeLSRippleCount(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.WAKEstate)));
    else
        wakeLSRippleRate(rec) = nan;
        wakeLSRippleCount(rec) = nan;
    end
    if ~isempty(ca1_ripples)
        wakeCA1RippleRate(rec) = length(Restrict(ca1_ripples.peaks,double(sleep.ints.WAKEstate)))./wakeTime;
    else
        wakeCA1RippleRate(rec) = nan;
    end
    if ~isempty(ca1_ripples) & ~isempty(ls_ripples)
        wakeCrossRegionCoupling(rec) = mean(coupling(find(InIntervals(ls_ripples.peaks,double(sleep.ints.WAKEstate)))));
    else
        wakeCrossRegionCoupling(rec) = nan;        
    end
    % get NREM rate
    if isfield(sleep.ints,'NREMstate')
        NREMTime = double(sum(diff(sleep.ints.NREMstate')));
        if NREMTime > 60
        if ~isempty(ls_ripples)
            NREMLSRippleRate(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.NREMstate)))./NREMTime;
            NREMLSRippleCount(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.NREMstate)));
        else
            NREMLSRippleRate(rec) = nan;
            NREMLSRippleCount(rec) = nan;
        end
        if ~isempty(ca1_ripples)
            NREMCA1RippleRate(rec) = length(Restrict(ca1_ripples.peaks,double(sleep.ints.NREMstate)))./NREMTime;
        else
            NREMCA1RippleRate(rec) = nan;
        end
        if ~isempty(ca1_ripples) & ~isempty(ls_ripples)
            NREMCrossRegionCoupling(rec) = mean(coupling(find(InIntervals(ls_ripples.peaks,double(sleep.ints.NREMstate)))));
        else
            NREMCrossRegionCoupling(rec) = nan;
        end
        else
            NREMLSRippleRate(rec) = nan;
            NREMCA1RippleRate(rec) = nan;
        end
    end
    % get REM rate
    if isfield(sleep.ints,'REMstate')
        REMTime = double(sum(diff(sleep.ints.REMstate')));
        if REMTime > 60
        if ~isempty(ls_ripples)
            REMLSRippleRate(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.REMstate)))./REMTime;
            REMLSRippleCount(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.REMstate)));
        else
            REMLSRippleRate(rec) = nan;
            REMLSRippleCount(rec) = nan;
        end
        if ~isempty(ca1_ripples)
            REMCA1RippleRate(rec) = length(Restrict(ca1_ripples.peaks,double(sleep.ints.REMstate)))./REMTime;
        else
            REMCA1RippleRate(rec) = nan;
        end
        if ~isempty(ca1_ripples) & ~isempty(ls_ripples)
            REMCrossRegionCoupling(rec) = mean(coupling(find(InIntervals(ls_ripples.peaks,double(sleep.ints.REMstate)))));
        else
            REMCrossRegionCoupling(rec) = nan;
        end
        end
    end
%     end
    
   cd /home/david/datasets/ripples_LS 
% cd E:\datasets\ripples_LS
   clear coupling
end


%% plotting
subplot(2,2,1)
errorbar(1,nanmean(wakeLSRippleRate(idx)),nanstd(wakeLSRippleRate(idx)),'.r')
hold on
errorbar(2,nanmean(NREMLSRippleRate(idx)),nanstd(NREMLSRippleRate(idx)),'.r')
errorbar(3,nanmean(REMLSRippleRate(idx)),nanstd(REMLSRippleRate(idx)),'.r')

errorbar(1.2,nanmean(wakeCA1RippleRate),nanstd(wakeCA1RippleRate),'.k')
hold on
errorbar(2.2,nanmean(NREMCA1RippleRate),nanstd(NREMCA1RippleRate),'.k')
errorbar(3.2,nanmean(REMCA1RippleRate),nanstd(REMCA1RippleRate),'.k')
axis([0 4 0.0001 .2])
set(gca,'yscale','log')
subplot(2,2,2)
errorbar(1,nanmean(wakeCrossRegionCoupling(idx)),nanstd(wakeCrossRegionCoupling(idx)),'.r')
hold on
errorbar(2,nanmean(NREMCrossRegionCoupling(idx)),nanstd(NREMCrossRegionCoupling(idx)),'.k')
errorbar(3,nanmean(REMCrossRegionCoupling(idx)),nanstd(REMCrossRegionCoupling(idx)),'.b')
