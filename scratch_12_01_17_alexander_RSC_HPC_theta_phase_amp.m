load('thetaPeakSpectDataForDavid.mat')

[b a] = butter(3,[4/626 12/625],'bandpass');
theta = FiltFiltM(b,a,hLFP);
[pks locs] = findpeaks(theta);
r = randperm(length(locs));
hspecs = zeros(200,1001);
rspecs = zeros(200,1001);
locs(1:4) = [];  % cut out the first few cycles
locs(end-4:end) = []; % cut out the last few cycles

for i=1:length(locs)
start = locs(r(i))-500;
stop = locs(r(i))+500;
[freqs,t,spec] = WaveSpec(hLFP(start:stop),[1 200],200,3,1/1000,'lin');
hspecs = hspecs + abs(spec);
[freqs,t,spec] = WaveSpec(rLFP(start:stop),[1 200],200,3,1/1000,'lin');
rspecs = rspecs+ abs(spec);
subplot(2,2,1)
imagesc(squeeze((hspecs)));
subplot(2,2,2)
imagesc(squeeze((rspecs)));
pause(.01)
end
for i=1:200
    rspecs_zscored(i,:) = zscore(rspecs(i,:));
    hspecs_zscored(i,:) = zscore(hspecs(i,:));
end
subplot(2,2,3)
imagesc(hspecs)
subplot(2,2,4)
imagesc(rspecs)
