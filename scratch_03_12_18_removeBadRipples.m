d = dir('*cgm*');
for i=1:length(d)
    cd(d(i).name)
    ripples = bz_LoadEvents(pwd,'CA1Ripples');
    scatter(ripples.data.peakAmplitude,ripples.data.peakFrequency,'.')
    if ~isempty(ripples)
    g = ginput();
    bad = find(ripples.data.peakAmplitude>g(1))

    ripples.timestamps(bad,:) = [];
    ripples.peaks(bad,:) = [];
    ripples.peakNormedPower(bad,:) = [];
    ripples.maps.ripples(bad,:) = [];
    ripples.maps.frequency(bad,:) = [];
    ripples.maps.phase(bad,:) = [];
    ripples.maps.amplitude(bad,:) = [];
    ripples.data.peakFrequency(bad,:) = [];
    ripples.data.peakAmplitude(bad,:) = [];
    ripples.data.duration(bad,:) = [];

    ripples
    pause

    save([d(i).name '.CA1Ripples.events.mat'],'ripples')
    end

    cd ..
end