for ff = 1:5
    
    isa = fillmissing((data{ff}.isig_levels),'nearest');
    isa = highpass(isa,dayHz,nyquist*2,'steepness',.5);

    minuteHz = 1/(60);
    hourHz = 1/(60*60);
    dayHz = 1/(24*60*60);

    lfp.samplingRate = minuteHz/5;
    lfp.data = isa';
    lfp.timestamps = data{ff}.absTime';
    [spec] = bz_WaveSpec(lfp,'frange',[dayHz minuteHz/10],'ncyc',3,'space','lin','nfreqs',1000);

    spec.data(isnan(data{ff}.isig_levels),:) = nan;

    plot(spec.freqs,nanmean(abs(spec.data)).*spec.freqs,'k')
    line([minuteHz minuteHz]/10,[.0001 .0006])
    line([hourHz hourHz],[.0001 .0006])
    line([hourHz hourHz]*2,[.0001 .0006],'color','r') % 30 min
    
end