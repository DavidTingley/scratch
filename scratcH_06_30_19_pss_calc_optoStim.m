% fileinfo = dir('analogin.dat');
% num_samples = fileinfo.bytes/(1*2);
% fid = fopen('analogin.dat','r');
% v = fread(fid,[1,num_samples],'uint16');
% fclose(fid)
% v = v*.000050354;

% analogin = downsample(v(1,:),16);
% 
%% specific for A32 and A45?
lfp.data = lfp.data(:,1);
lfp.channels = lfp.channels(1);
%%
% 
% [pks pulseIdx] = findpeaks(fastrms(analogin,20),'minPeakheight',.5);
pulseTimes = lfp.timestamps(pulseIdx);
[b a] = butter(3,50/625,'high');

for p = 1:length(pulseIdx)
    stims(p,:) = fastrms(analogin(pulseIdx(p)-625:pulseIdx(p)+625),10);
    responses(p,:) = lfp.data(pulseIdx(p)-625:pulseIdx(p)+625);
    responses_shuf(p,:) = lfp.data(-2*1250-625 + [pulseIdx(p)-625:pulseIdx(p)+625]);
    r(p,:) = filtfilt(b,a,double(responses(p,:)));
    r_shuf(p,:) = filtfilt(b,a,double(responses_shuf(p,:)));
end

for i=1:length(pulseIdx)
[spec]= bz_WaveSpec(responses(i,:)','samplingRate',1250,'space','lin','frange',[1 300]);
s(i,:,:) = abs(spec.data);
end
whos s
ss = (squeeze(mean(s))')
whos s ss
imagesc(ss)
caxis
for i=1:100
ss(i,:) = zscore(ss(i,:));
end
imagesc(ss)


[specslope] = bz_PowerSpectrumSlope(lfp,4,1);


for ii=1:45
    idx = find(clu==ii);
%     imagesc(r(idx,:))
%     caxis([-50 50])
%     [x y] = ginput();
%     for i=1:length(idx)
%     pow(i) = max(abs(r(idx(i),x(1):x(2))));
%     end
    stimDur(ii) = sum(find(nanmean(stims(idx,:))>.05));
    nStims(ii) = length(idx);
    for i=1:length(idx)
    pow(i) = max(smooth(abs(r(idx(i),625-100:625+100)),40));
    end
    for i=1:length(idx)
    p = pulseTimes(idx(i));
    [a b] = min(abs(specslope.timestamps-p));
    pss(i) = specslope.data(b,1);
    end

    pow(pow==0)=nan;
    power{ii} = (pow);
    powspec{ii} = (pss);
    [pss_pow_corr(ii) pv(ii)]= corr(pow',pss','rows','complete'); clear pow pss

end



for p = 1:length(pulseIdx)
[stimAmp(p) b] = max(stims(p,:));
end

idx = 1:length(stimAmp); %find(stimAmp>.25);
clear psd*
parfor p = 1:length(idx)
h = spectrum.periodogram;
hpsd = psd(h,double(responses(idx(p),125:1000)),'Fs',1250);
ps = zeros(length(0:.125:250),1);
c=1;
for hz = 0:.125:250
t = hz-.125 <= hpsd.Frequencies & hpsd.Frequencies<hz;
ps(c) = mean(log(hpsd.Data(t)));
c=c+1;
end
psd_actual(p,:) = ps;
end
parfor p = 1:length(idx)
h = spectrum.periodogram;
hpsd = psd(h,double(responses_shuf(idx(p),125:1000)),'Fs',1250);
ps = zeros(length(0:.125:250),1);
c=1;
for hz = 0:.125:250
t = hz-.125 <= hpsd.Frequencies & hpsd.Frequencies<hz;
ps(c) = mean(log(hpsd.Data(t)));
c=c+1;
end
psd_shuf(p,:) = ps;
end
id = find(~isnan(psd_actual(1,:)));
idx = find(stimAmp>.25);
boundedline(vec(id),nanmean(psd_shuf(idx,id)),nanstd(psd_shuf(idx,id)),'r')
hold on
boundedline(vec(id),nanmean(psd_actual(idx,id)),nanstd(psd_actual(idx,id)),'k')

for i=1:2001
[a(i) b(i)] = ttest2(psd_actual(idx,i),psd_shuf(idx,i));
end


stimThresh = [.25 .35 .45 .55 .65];
for s = 1:5
    idx = find(stimAmp>stimThresh(s));
    for i=1:2001
        [a(s,i) b(s,i)] = ttest2(psd_actual(idx,i),psd_shuf(idx,i));
    end
end
whos a b
plot(b(:,id)')







