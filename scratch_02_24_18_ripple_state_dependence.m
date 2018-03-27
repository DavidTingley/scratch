d = dir('*201*');
count=1;
for rec = 1:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    sessionInfo.FileName
    ls_ripples = bz_LoadEvents(pwd,'LSRipples');
    ca1_ripples = bz_LoadEvents(pwd,'CA1Ripples');
    sleep = bz_LoadStates(pwd,'SleepState');
    
    % get wake rate
    wakeTime = double(sum(diff(sleep.ints.WAKEstate')));
    if ~isempty(ls_ripples)
        wakeLSRippleRate(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.WAKEstate)))./wakeTime;
    else
        wakeLSRippleRate(rec) = nan;
    end
    if ~isempty(ca1_ripples)
        wakeCA1RippleRate(rec) = length(Restrict(ca1_ripples.peaks,double(sleep.ints.WAKEstate)))./wakeTime;
    else
        wakeCA1RippleRate(rec) = nan;
    end
    
    % get NREM rate
    NREMTime = double(sum(diff(sleep.ints.WAKEstate')));
    if NREMTime > 600
    if ~isempty(ls_ripples)
        NREMLSRippleRate(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.NREMstate)))./NREMTime;
    else
        NREMLSRippleRate(rec) = nan;
    end
    if ~isempty(ca1_ripples)
        NREMCA1RippleRate(rec) = length(Restrict(ca1_ripples.peaks,double(sleep.ints.NREMstate)))./NREMTime;
    else
        NREMCA1RippleRate(rec) = nan;
    end
    else
        NREMLSRippleRate(rec) = nan;
        NREMCA1RippleRate(rec) = nan;
    end
    % get REM rate
    REMTime = double(sum(diff(sleep.ints.WAKEstate')));
    if ~isempty(ls_ripples)
        REMLSRippleRate(rec) = length(Restrict(ls_ripples.peaks,double(sleep.ints.REMstate)))./REMTime;
    else
        REMLSRippleRate(rec) = nan;
    end
    if ~isempty(ca1_ripples)
        REMCA1RippleRate(rec) = length(Restrict(ca1_ripples.peaks,double(sleep.ints.REMstate)))./REMTime;
    else
        REMCA1RippleRate(rec) = nan;
    end
   cd /home/david/datasets/ripples_LS 
end


%% plotting

errorbar(1,nanmean(wakeLSRippleRate),nanstd(wakeLSRippleRate)./sqrt(93)*3,'.r')
hold on
errorbar(2,nanmean(NREMLSRippleRate),nanstd(NREMLSRippleRate)./sqrt(93)*3,'.r')
errorbar(3,nanmean(REMLSRippleRate),nanstd(REMLSRippleRate)./sqrt(93)*3,'.r')

errorbar(1.2,nanmean(wakeCA1RippleRate),nanstd(wakeCA1RippleRate)./sqrt(93)*3,'.k')
hold on
errorbar(2.2,nanmean(NREMCA1RippleRate),nanstd(NREMCA1RippleRate)./sqrt(93)*3,'.k')
errorbar(3.2,nanmean(REMCA1RippleRate),nanstd(REMCA1RippleRate)./sqrt(93)*3,'.k')
axis([0 4 0.0001 .2])
set(gca,'yscale','log')

