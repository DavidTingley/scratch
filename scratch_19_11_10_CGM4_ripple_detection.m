
d = dir('*CGM*');
ripChans = [];

for f = 1:length(d) % left at 5
    cd(d(f).name)
    file = dir('*CA1Ripples*');
    orig = load(file.name);
    plot(f,orig.ripples.stdev,'.k')
    hold on
    lfp = bz_GetLFP([16 60]); % 26, or 8 or 6, for 288 
    lfp.data(:,1) = orig.ripples.detectorinfo.detectionparms.lfp;

    [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,'noise',lfp.data(:,2),'stdev',30000,'emgthresh',0,'thresholds',[2 4],'durations',[15 250]);
    ripples.stdev
%     [b a] = butter(4,[120/625 200/625],'bandpass');
%     filt  = filtfilt(b,a,double(lfp.data(:,1)));
% 
%     [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(filt,lfp.timestamps,ripples);
% 
%     save([d(f).name '.ripples.events.mat'],'ripples')
    cd ..
end


