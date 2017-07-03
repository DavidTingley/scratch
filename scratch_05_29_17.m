d = dir('*')
for i=3:length(d)
    if d(i).isdir
        cd(d(i).name)
            r = dir('*ripplesALL*');
            if ~isempty(r)
                load(r.name)
                if isfield(ripples,'detectorparms')
                    ripples.detecorParams = ripples.detectorparms;
                ripples = rmfield(ripples,'detectorparms');
                end
                xml = LoadParameters;
                rippleChan = ripples.rippleChan;
            refChan = ripples.refChan;
           stdev = (ripples.stdev)*1.5;
 csv = readtable('/mnt/packrat/userdirs/david/zpool2/metadata.csv');
        temp = pwd;
        temp = strsplit(temp,'/');
        animal = temp{length(temp)-1};
        recording = temp{length(temp)};
        temp{length(temp)} = []; temp = removeEmptyCells(temp);
        t = strsplit(recording,'_');
%         t=str2num(t{end});
%         if t == 10 | t==9 | t == 7 | t == 6 | t == 3
%             rippleChan = 12;
%         elseif t == 8
%             rippleChan = 10    
%         elseif t == 6 | t == 4 | t == 5
%             rippleChan = 4
%         end
         lfp = bz_GetLFP([rippleChan refChan]);
        animalFolder = strjoin(temp,'/');
            [b a] = butter(4,[120/(xml.lfpSampleRate./2) 200/(xml.lfpSampleRate./2)],'bandpass');
            for i=1:2
                lfp.filt(:,i) = FiltFiltM(b,a,double(lfp.data(:,i)));
            end

            % need to get intervals of SWS here...

                load([xml.FileName '.SleepState.states.mat'])

                disp(['finding ripples with ch #: ' num2str(rippleChan) ', and noise ch #: ' num2str(refChan)])
                [ripples] = bz_FindRipples(lfp.data(:,1),lfp.timestamps,'noise',...
                    lfp.data(:,2),'durations',[20 300] ,'saveMat',false,'frequency',lfp.samplingRate);
                [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripples,'frequency',lfp.samplingRate);

                ripples.rippleChan = rippleChan;
                ripples.refChan = refChan;

                save(['/' animalFolder '/' recording '/' xml.FileName '.ripplesALL.event.mat'],'ripples','-v7.3')

                %% parse by NREM
                ind = [];
                for s = 1:size(SleepState.ints.NREMstate,1)
                    temp = FindInInterval(ripples.peaks,double(SleepState.ints.NREMstate(s,:)));
                    if ~isempty(temp)
                    ind = [ind,temp(1):temp(2)];
                    end
                end
                ripplesNREM.times = ripples.times(ind,:);
                ripplesNREM.peaks = ripples.peaks(ind);
                ripplesNREM.stdev = ripples.stdev;
                ripplesNREM.noise = ripples.noise;
                ripplesNREM.peakNormedPower = ripples.peakNormedPower(ind);
                ripplesNREM.detectorParams = ripples.detectorParams;
                ripplesNREM.detectorName = ripples.detectorName;
                if length(ripplesNREM.times)>2
                    [ripplesNREM.maps,ripplesNREM.data,ripplesNREM.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripplesNREM,'frequency',lfp.samplingRate);
                    save(['/' animalFolder '/' recording '/' xml.FileName '.ripplesNREM.event.mat'],'ripplesNREM','-v7.3')
                end
                %% parse by WAKE
                ind = [];
                for s = 1:size(SleepState.ints.WAKEstate,1)
                    temp = FindInInterval(ripples.peaks,double(SleepState.ints.WAKEstate(s,:)));
                    if ~isempty(temp)
                    ind = [ind,temp(1):temp(2)];
                    end
                end
                ripplesWAKE.times = ripples.times(ind,:);
                ripplesWAKE.peaks = ripples.peaks(ind);
                ripplesWAKE.stdev = ripples.stdev;
                ripplesWAKE.noise = ripples.noise;
                ripplesWAKE.peakNormedPower = ripples.peakNormedPower(ind);
                ripplesWAKE.detectorParams = ripples.detectorParams;
                ripplesWAKE.detectorName = ripples.detectorName;
                if length(ripplesWAKE.times)>2
                    [ripplesWAKE.maps,ripplesWAKE.data,ripplesWAKE.stats] = bz_RippleStats(lfp.filt(:,1),lfp.timestamps,ripplesWAKE,'frequency',lfp.samplingRate);
                    save(['/' animalFolder '/' recording '/' xml.FileName '.ripplesWAKE.event.mat'],'ripplesWAKE','-v7.3')
                end
            end
                
        cd ..
    end
end
