d = dir('*');

for i=5:-1:3
    cd(d(i).name)
    
    dd = dir('*_21*');
    for j=1:length(dd)
       cd(dd(j).name)
       dors = bz_LoadEvents(pwd,'DorsalRipples');
       vent = bz_LoadEvents(pwd,'VentralRipples');
       if ~isempty(dors) & ~isempty(vent)
           bad = [];

           for r = 1:length(dors.peaks)
              [a b] = min(abs(dors.peaks(r)-vent.peaks));
              if a > 0.06
                  bad = [bad, r];
              end
           end

           ripples = dors;

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
            clear bad

           save([dd(j).name '.GlobalRipples.events.mat'],'ripples') 
       end
       cd ..
    end   
    cd ..
    d(i).name
end
