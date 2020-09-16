cd('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data')
d = dir('_*');

for a = 1:8
    refract{a} = [];
    rips{a} = [];
    freqs{a} = [];
    amps{a} = [];
    zamps{a} = [];
    durs{a} = [];
    ripRate{a} = [];
    thetaR{a} = [];
    
    cd('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data')
    load(d(a).name,'count','isig_levels','rec','idx')
    animal = strsplit(d(a).name,'_');
    animal = strsplit(animal{2},'.');
    animal = animal{1};
    cd(['D:\datasets\glucose\' animal])
    recs = dir(['*' animal '_*']);
    
    for r = 1:length(recs)
        cd(recs(r).name)
        f = find(rec(idx)==r);
        sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
        if length(f) > 10 & exist([sessionInfo.FileName '.CA1Ripples.events.mat'])
            load([sessionInfo.FileName '.CA1Ripples.events.mat'])
            load([sessionInfo.FileName '.thetaResid.mat'])
            
            refract{a} = [refract{a};diff(ripples.peaks)];
            temp = zeros(length(theta_resid)*1250,1);
            temp(round(ripples.peaks*1250)) = 1;
%             ripRate{a} = [ripRate{a};makelength(fastrms(temp,2500),length(theta_resid))'];clear temp
            ripRate{a} = [ripRate{a};makeLength(fastrms(temp,125),round(length(temp)./50))'];clear temp
            rips{a} = [rips{a};ripples.maps.ripples];
            freqs{a} = [freqs{a};ripples.data.peakFrequency];
            amps{a} = [amps{a};(ripples.data.peakAmplitude)];
            zamps{a} = [zamps{a};nanZscore(ripples.data.peakAmplitude)];
            durs{a} = [durs{a};ripples.data.duration];
            thetaR{a} = [thetaR{a};theta_resid'];
        end  
        cd ..
    end
end


[pxx_filt f] = periodogram(ripRate{a},[],256,1);

