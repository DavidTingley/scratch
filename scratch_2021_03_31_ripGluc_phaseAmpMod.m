clear all
orig_files = dir('_*mat');
new_files = dir('revision/*mat');

% cgm_files = dir('stim/cgm*mat');
cgm_files = [orig_files;new_files];

for j=1:length(cgm_files)
    if ~strcmp(cgm_files(j).name,'CGM36.mat')
        load([cgm_files(j).folder '/' cgm_files(j).name],'absTime','nyquist','*Hz','count','stimRate','isig_levels','isa','idx','spSlope','emgSig','mov','theta_z','states')
%         count = stimRate; 
        if ~strcmp(cgm_files(j).name,'CGM47.mat') & ~strcmp(cgm_files(j).name,'CGM50.mat')
            isig = nanZscore([0, diff(isig_levels)]);
        else
            isig = nanZscore(isa);
        end
        if strcmp(cgm_files(j).name,'ros.mat')
            isig = nanZscore(isa);
        end
        
        remove = find(~ismember(1:length(count),idx));

        
        lfp.samplingRate = 1/(60*5);
        
        [b a] = butter(4,[ dayHz / nyquist],'high');
        
        count(remove)=nan;
        count(count>250)=nan;
%         count(count==0)=nan;
        isig(remove)=nan;

%         isig = fillmissing(isig,'nearest');
%         count = fillmissing(count,'nearest');
        cut = find(isnan(count) | isnan(isig));
        absTime(cut)=[];
        count(cut) =[];
        isig(cut) = [];
        lfp.timestamps = absTime';
        lfp.data(:,1) = count; %filtfilt(b,a,count);
        lfp.refCh = circshift(count',-2);
        lfp.data(:,2) = isig; %filtfilt(b,a,isig);
        
        lfp.channels = [1 2];
        %% phase/amplitude
        freqRange = dayHz*2:dayHz/2:nyquist/1.5;
        try
            [comodulogram_ripGluc] = bz_CFCPhaseAmp(lfp,freqRange,freqRange,'phaseCh',1,'ampCh',2,'makePlot',false);
%             title('ripGluc')
            [comodulogram_glucRip] = bz_CFCPhaseAmp(lfp,freqRange,freqRange,'phaseCh',2,'ampCh',1,'makePlot',false);

%             title('glucRip')
            ripGluc(j,:,:)=(comodulogram_ripGluc.comod')';
            glucRip(j,:,:)=(comodulogram_glucRip.comod')';
        catch
            ripGluc(j,:,:)=nan(length(freqRange)-1,length(freqRange)-1);
            glucRip(j,:,:)=nan(length(freqRange)-1,length(freqRange)-1);
        end
        %% power/power
        parms.frange = [dayHz nyquist];
        parms.nfreqs = 200;
        [comod] = bz_Comodulogram(lfp,parms);
        ripGluc_powpow(j,:,:) = comod.corrs{2};
        
        %% plotting
        [a b] = min(abs(freqRange-hourHz));
        [a bb] = min(abs(freqRange-hourHz*2));
        subplot(3,2,1)
        imagesc((squeeze(nanmean(ripGluc))')')
        xline(b);
        xline(bb);
%         xline(84);
        title('ripGluc')
        subplot(3,2,2)
        imagesc((squeeze(nanmean(glucRip))')')
        xline(b);
        xline(bb);
%         xline(84);
        title('glucRip')
        subplot(3,2,3)
        imagesc((squeeze((ripGluc(j,:,:)))))
        xline(b);
        xline(bb);
%         xline(84);
        title('ripGluc')
        subplot(3,2,4)
        imagesc((squeeze((glucRip(j,:,:)))))
        xline(b);
        xline(bb);
        title('glucRip')
        subplot(3,2,5)
        plot(lfp.data(:,1))
        subplot(3,2,6)
        
        
        plot(lfp.data(:,2))
        title(j)
%         
        pause(1)
        
        ang = angle(hilbert(isig));
        c=1;
        for t = -pi:.1:pi
           id = find(ang>t-.1 & ang<t+.1);
           if ~isempty(id)
               phaseLocked(j,c) = nanmedian(count(id));
           else
               phaseLocked(j,c)=0;
           end
           c=1+c;
        end
        c=1;
        for i=1:length(count)
            for jj=1:count(i)
                a(c) = ang(i);
                c=1+c;
            end
        end
        circ_locking(j) = circ_rtest(a);clear a
        
        clear lfp remove isig
    end
end


imagesc(squeeze(nanmedian(ripGluc_powpow)))
caxis([-.3 .3])
hold on
plot([0 200],[0 200],'k')

xlabel('SPW-R rate frequency')
ylabel('glucose frequency')
title('power/power correlation (white-6hours; red-60min; black-15min)')
colorbar
xline(56,'w')
xline(127,'r')
xline(182,'k')