cd('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data')
d = dir('_*');

for a = 1:8
    cd('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data')
    load(d(a).name,'count','isig_levels','rec','idx')
    animal = strsplit(d(a).name,'_');
    animal = strsplit(animal{2},'.');
    animal = animal{1};
    cd(['D:\datasets\glucose\' animal])
    recs = dir(['*' animal '_*']);
    
    for r = 1:length(recs)
        cd(recs(r).name)
        f = find(rec==r);
        sessionInfo = bz_getSessionInfo(pwd,'noprompts',true);
        
        if length(f) > 10 & exist([sessionInfo.FileName '.CA1Ripples.events.mat'])
            load([sessionInfo.FileName '.CA1Ripples.events.mat'])
%             lfp = bz_GetLFP(1);
            lfp.data = ripples.detectorinfo.detectionparms.lfp;
            lfp.samplingRate = 1250;
            lfp.timestamps = [0:1/1250:(length(lfp.data)-1)./1250]';
            clear ripples;
            spec = bz_WaveSpec(lfp,'frange',[1 200],'nfreqs',20,'ncyc',3,'space','lin','downsampleout',1250);        
            clear lfp
            
            glucFlux = [0,diff(isig_levels(f))];
            gluc = [(isig_levels(f))];
            glucFlux_int=makelength(glucFlux,length(spec.data));
            gluc_int=makelength(gluc,length(spec.data));

            c=1; 
            for j=-3600:60:3600
                for i=1:20
                cc(r,i,c) = corr(circshift(abs(spec.data(:,i)),j),glucFlux_int','rows','complete');
                cc_gluc(r,i,c) = corr(circshift(abs(spec.data(:,i)),j),gluc_int','rows','complete');
                end
                c=c+1;
                subplot(2,2,1)
                imagesc(zscore(squeeze(nanmean(cc,1))',[],1)');
                subplot(2,2,2)
                imagesc(zscore(squeeze(nanmean(cc_gluc,1))',[],1)');
                
                pause(.1)
            end  
        else
            cc(r,1:20,1:121)=nan;
            cc_gluc(r,1:20,1:121)=nan;
        end
        r
        cd ..
    end
    corrFlux{a}=cc; clear cc
    corrGluc{a}=cc_gluc; clear cc_gluc
    a
end