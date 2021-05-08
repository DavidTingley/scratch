d = dir('*STZ*');
for i = 1:length(d)
    if d(i).isdir
        cd(d(i).name)
        ripples = bz_LoadEvents(pwd,'CA1Ripples');
        if ~isempty(ripples)
            subplot(2,1,1)
        scatter(ripples.data.peakAmplitude,ripples.data.peakFrequency,'.')
        subplot(2,1,2)
        scatter(ripples.data.peakAmplitude,ripples.data.duration,'.')
        g = ginput();
        bad = find(ripples.data.peakAmplitude>g(1));

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
        if length(bad)>1
            save([d(i).name '.CA1Ripples.events.mat'],'ripples')
        end
        else
            disp([d(i).name ' recording has no ripples?'])
        end

        cd ..
    end
end