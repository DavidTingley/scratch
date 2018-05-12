d = dir('*201*');

for rec=1:length(d)
    cd(d(rec).name)
    if exist([d(rec).name '.behavior.mat'])
        load([d(rec).name '.behavior.mat'])
        ls.ripples = bz_LoadEvents(pwd,'LSRipples');
        ca1.ripples = bz_LoadEvents(pwd,'CA1Ripples');
        if ~isempty(ls.ripples)
            lsRate = zeros(length(behavior.events.trialIntervals),6000);
        else
            lsRate = nan(length(behavior.events.trialIntervals),6000);
        end
        if ~isempty(ls.ripples)
            ca1Rate = zeros(length(behavior.events.trialIntervals),6000);
        else
            ca1Rate = nan(length(behavior.events.trialIntervals),6000);
        end
        
        for trial = 1:length(behavior.events.trialIntervals)
            ts = behavior.events.trialIntervals(trial,2);
            if ~isempty(ls.ripples)
                ind1 = find(abs(ls.ripples.peaks-ts)<3);
                lsRate(trial,round((ls.ripples.peaks(ind1)-ts+3)*1000)) = 1;
            end
            if ~isempty(ca1.ripples)
                ind2 = find(abs(ca1.ripples.peaks-ts)<3);
                ca1Rate(trial,round((ca1.ripples.peaks(ind2)-ts+3)*1000)) = 1;
            end           
        end       
        
        lsRippleRate(rec,:) = mean(lsRate);
        ca1RippleRate(rec,:) = mean(ca1Rate);
        numTrials(rec) = trial;
        subplot(3,2,1)
        imagesc(lsRippleRate);
        subplot(3,2,2)
        imagesc(ca1RippleRate);
        
        subplot(3,2,3)
        plot(nanmean(lsRippleRate))
        subplot(3,2,5)
        plot(smooth(nanmean(lsRippleRate),250))
        subplot(3,2,4)
        plot(nanmean(ca1RippleRate))
        subplot(3,2,6)
        plot(smooth(nanmean(ca1RippleRate),250))
        pause(.01)
        clear behavior lsRate ca1Rate ls ca1 ind1 ind2 ts
    end
    rec
    cd /home/david/datasets/ripples_LS/
end