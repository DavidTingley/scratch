d = dir('*201*')
for rec = 5:length(d)
    cd(d(rec).name);
load([d(rec).name '.jump.behavior.mat'])
sessionInfo = bz_getSessionInfo;
lfp = bz_GetLFP(sessionInfo.ca1);
if isfield(sessionInfo,'ls') 
septal = bz_GetLFP(sessionInfo.ls);
end
if ~isempty(sessionInfo.ca1)
    ints = mean(behavior.events.trialIntervals,2);
    for i=1:length(behavior.events.trialConditions)
    if ~isnan(ints(i))
    %ca1
    [aa bb] = min(abs(lfp.timestamps-ints(i)));
    [freqs,t,spec] = WaveSpec(double(lfp.data(bb-1250*3:bb+1250*3)),[1 200],200,3,1/1250,'lin');
    specs(i,:,:) = abs(spec);
    subplot(2,2,1)
    imagesc(squeeze(mean(specs)))
    subplot(2,2,2)
    for j=1:200
    ss(i,j,:) = zscore(abs(specs(i,j,:)));
    end
    end
    end
    s = squeeze(median(ss));
    for i=1:200
    sss(i,:) = smooth(zscore(s(i,:)),50);
    end
    lfp_wavelet{rec} = sss;
    imagesc(squeeze(mean(ss)))
    %ls
    if isfield(sessionInfo,'ls') 
        if ~isempty(sessionInfo.ls)
    for i=1:length(behavior.events.trialConditions)
    if ~isnan(ints(i))
    [freqs,t,spec] = WaveSpec(double(septal.data(bb-1250*3:bb+1250*3)),[1 200],200,3,1/1250,'lin');
    specs(i,:,:) = abs(spec);
    subplot(2,2,3)
    imagesc(squeeze(mean(specs)))
    subplot(2,2,4)
    for j=1:200
    ss(i,j,:) = zscore(abs(specs(i,j,:)));
    end
    end
    end
    imagesc(squeeze(mean(ss)))
    s = squeeze(median(ss));
    for i=1:200
    sss(i,:) = smooth(zscore(s(i,:)),50);
    end
    septal_wavelet{rec} = sss;
        end
    end
    pause(.01)

    clear s ss sss spec specs
end
    
    cd /home/david/datasets/jumpDataset/
end