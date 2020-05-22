clear all
cd /home/david/Dropbox/Documents/pubs/inProgress/glucose

% parpool(4)
animal = 'cgm8'; 
addpath(genpath(pwd));
read_CGM_data
clear a b dat nSteps timeStep ts ISIG_vals datatype dates times % clean up
read_stimLog_data
cd(['/mnt/packrat/userdirs/david/zpool4/' animal])
% cd /mnt/nyuShare/Buzsakilabspace/david/zpool2/CGM2
d = dir(['*' animal '*']);

for i=1:length(d)
    if d(i).isdir
        cd(d(i).name);
        ripples{i} = bz_LoadEvents(pwd,'CA1Ripples');
        if ~isempty(ripples{i})
        ripples{i}.data.zScoredAmps = zscore(ripples{i}.data.peakAmplitude);
        end
        if ~isempty(ripples{i})
            %% load some stuff
            emgFromLFP{i} = load([d(i).name '.EMGFromLFP.LFP.mat']);
            emgFromLFP{i}.EMGFromLFP.data = single(emgFromLFP{i}.EMGFromLFP.data);
            
            [zt_rec{i} rt{i} at{i}] = bz_getZeitgeberTime(d(i).name,20000,20000);
            SleepState{i} = bz_LoadStates(pwd,'SleepState');

            %% add accelerometer data
            if ~exist([d(i).name '.movement.mat'])
                num_channels = 3;
                fid = fopen('auxiliary.dat', 'r');
                f = dir('auxiliary.dat');
                v = 1;
                chunk = 1;
                while ~isempty(v)
                    v = single(fread(fid, [num_channels, 100000], 'uint16=>uint16',6*5));
                    v = single(v * 0.0000374); % convert to volts
    %                 chunk
                    if ~isempty(v)
                    for t = 1:size(v,1) 
                        diffs(t,:)= single(diff([v(t,:)]')); % smooth and interp to seconds
                    end                
                    % exclude bad aux channel...
                    if any(std(diffs')>10*median(std(diffs')))
                        diffs(find(std(diffs')>10*median(std(diffs'))),:) = nan;
                    end
                    movement{i}(chunk:chunk+length(v)-2) = nanmean(abs(diffs)); 
                    if chunk ~= 1 
                    movement{i}(chunk+length(v)-1) = nanmean(abs(lastV-v(:,1)));
                    end
                    lastV = v(:,end);
                    clear diffs 
                    end
                    chunk = chunk + 100000;
    %                 movement{i} = sum(abs(diffs)); clear diffs   
                end
                movement_s = makeLength(movement{i},f.bytes/6/5000);
                clear v diffs lastV
                fclose(fid);
                save([d(i).name '.movement.mat'],'movement_s')
                movement{i} = movement_s; clear movement_s;
            else
                load([d(i).name '.movement.mat'],'movement_s')
                movement{i} = movement_s; clear movement_s;
            end
            
            %% get stim times
            if ~exist([d(i).name '.stimTimes.mat'])
%                 if exist('digitalin.dat')
%                     fileinfo = dir('digitalin.dat');
%                     num_samples = fileinfo.bytes/2;
%                     fid = fopen('digitalin.dat','r');
%                     dig = fread(fid,num_samples,'uint16=>uint16');
%                     digSort = bitand(dig,2^1)>0; clear dig;
%                     for chunk = 1:100000:length(digSort)
%                         if chunk + 100000 < length(digSort)
%                             di(chunk:chunk+99999) = int16(diff(digSort(chunk:chunk+100000)));
%                             di(end+1) = digSort(chunk+100000)-digSort(chunk+100001);
%                         else
%                             di(chunk:length(digSort)-1) = int16(diff(digSort(chunk:end)));
%                             di(end+1) = digSort(length(digSort)-1)-digSort(end);
%                         end
%                     end
%                     clear digSort
%                     pos = find(di==1); clear di dig 
%             %         neg = find(d==-1);
%                     stimTimes_s = pos/20000;  % relative timing w/ in recording
%                     clear pos
%                 else
%                    stimTimes_s = []; 
%                 end
%                 save([d(i).name '.stimTimes.mat'],'stimTimes_s')
%                 stimTimes{i} = stimTimes_s; clear stimTimes_s
            else
                load([d(i).name '.responseMagnitudes.mat'])
                load([d(i).name '.stimTimes.mat'],'pul')
                stims_analog = pul{1}; clear pul
                
                inRec = find(at{i}(1) < stimTimes_logFile & at{i}(end) > stimTimes_logFile);
                c=1;
                stims = [];
                for t = 1:length(inRec)
                    if min(abs(stimTimes_logFile(inRec(t))-at{i})) < .005
                        stims(c)= stimTimes_logFile(inRec(t));
                        c=c+1;
                    end
                end
                stims = sort(stims);
                
                if length(stims)==length(stims_analog)
                    stimTimes{i} = stims_analog; clear stims stims_analog
                elseif length(stims) > length(stims_analog)
                    stimTimes{i} = stims; clear stims stims_analog
                elseif length(stims) < length(stims_analog)
                    stimTimes{i} = stims_analog; clear stims stims_analog
                else
                    error('wut')
                end
                responses{i} = responseMag;
            end
        
             %% add power spectrum
            if ~exist([d(i).name '.specSlope.mat'])
                lfp.data = ripples{i}.detectorinfo.detectionparms.lfp;
                lfp.timestamps = [0.00001:length(lfp.data)]'./1250;
                lfp.samplingRate = 1250;
                [specslope_singleRec,spec] = bz_PowerSpectrumSlope(lfp,20,10,'showfig',false,'nfreqs',10);
                save([d(i).name '.specSlope.mat'],'specslope_singleRec','spec')
                specslope{i} = specslope_singleRec; clear specslope_singleRec;
            else
                load([d(i).name '.specSlope.mat'])
                specslope{i} = specslope_singleRec; clear specslope_singleRec;
            end

            %% theta power
            if ~exist([d(i).name '.thetaPower.mat'])
                [bb aa] = butter(3,[6/625 11/625],'bandpass');
                temp = decimate(filtfilt(bb,aa,double(ripples{i}.detectorinfo.detectionparms.lfp)),1250);
                temp(1:10)=nan; temp(end-10:end)=nan; % kill those edges
                thetaPow_raw_s = fastrms(temp,60*5);
                thetaPow_z_s = nanZscore(thetaPow_raw_s);
                save([d(i).name '.thetaPower.mat'],'thetaPow*')
                thetaPow_raw{i} = thetaPow_raw_s; clear thetaPow_raw_s
                thetaPow_z{i} = thetaPow_z_s; clear thetaPow_z_s
            else
                load([d(i).name '.thetaPower.mat'],'thetaPow*_s')
                thetaPow_raw{i} = thetaPow_raw_s; clear thetaPow_raw_s
                thetaPow_z{i} = thetaPow_z_s; clear thetaPow_z_s
            end
                
        else
            emgFromLFP{i} = [];
            thetaPow_z{i} = [];
            thetaPow_raw{i} = [];
            specslope{i} = [];
            SleepState{i} = [];
            zt_rec{i} = [];
            rt{i} = [];
            at{i} = [];
        end
        cd(['/mnt/packrat/userdirs/david/zpool4/' animal])
%         cd /mnt/nyuShare/Buzsakilabspace/david/zpool2/CGM1
        
        % plot(cell2mat(at));
        % pause(.1)
    else
        ripples{i} = [];
        SleepState{i} = [];
        zt_rec{i} = [];
        rt{i} = [];
        at{i} = [];
    end
end


remove = [1:20];
% 1860:2419
isig_orig = isig_levels;

isig_levels(remove) = nan;

nyquist = 5./(24*60)/2;
dayHz = 1/(24*60*60);
hourHz = 1/(60*60);
minuteHz = 1/(60);
isa = fillmissing((isig_levels),'nearest');
% isa2 = highpass(isa,dayHz,nyquist*2);
[b a] = butter(4,[ dayHz / nyquist],'high');
isa = filtfilt(b,a,isa);
% isa = filtfilt(b,a,isa);
isa(isnan(isig_levels)) = nan;


[hi lo]=butter(4,[140/625 180/625],'bandpass');
c=1;
c_NREM = 1;
c_WAKE = 1;
c_stim = 1;
c_slope = 1;
stimAlignedCGM = [];
responseMag=[];
stimTime =[];
stimTimesAbs = [];
specSlope = [];
specSlopeTimes =[];
rippleAlignedCGM = [];
ripTime = [];
ripState = [];
offsets = [];
m = [];
zTime = [];
absolute = [];
amplitude = [];
duration = [];
frequency = [];
recording = [];
thetaPower_z = [];
thetaPower_raw = [];
state = [];
emg = [];

for ii=1:length(d)
    if  d(ii).isdir
if ~isempty(rt{ii}) & ~isempty(ripples{ii}) % & ~isempty(SleepState{ii})
    relTime = rt{ii}; 
    abTime = at{ii}; 
    rips = ripples{ii}; 
    sti = stimTimes{ii};
    if ~isempty(stimTimes{ii})
    stimAlignedCGM_t = nan(length(stimTimes{ii}),48);
%     responseMag_t = nan(length(stimTimes{ii}),1251);
    stimTime_t = nan(length(stimTimes{ii}),1);
    stimTimesAbs_t = nan(length(stimTimes{ii}),1);
    
    parfor s = 1:length(sti)
        [a b]= min(abs(relTime-sti(s)));
        [aa bb] = min(abs(abTime(b)-absTime));
        stimTimesAbs_t(s) = absTime(bb);
        
        offset = (abTime(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
        idx = find(absTime> abTime(b)-1.15741277113557e-05*60*60*2 & ...
            absTime < abTime(b)+1.15741277113557e-05*60*60*2);
        if length(idx)>48 & ~isempty(idx)
            idx = idx(1:48);
        elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
            idx(end+1) = idx(end) + 1;
        end
        if length(idx) == 48
            responseMag_t(s,:) = responses{ii}(s,:);
            stimAlignedCGM_t(s,:) = isa(idx);
            stimTime_t(s) = abTime(b);
        else
            responseMag_t(s,:) = responses{ii}(s,:);
            stimAlignedCGM_t(s,:) = nan(48,1);
            stimTime_t(s) = abTime(b);
        end
    end
    stimAlignedCGM = [stimAlignedCGM;stimAlignedCGM_t];
    responseMag = [responseMag;responseMag_t];
    stimTime = [stimTime;stimTime_t];
    stimTimesAbs = [stimTimesAbs;stimTimesAbs_t];
    end
    
    specSlope_t = nan(length(specslope{ii}.timestamps),1);
    specSlopeTimes_t = nan(length(specslope{ii}.timestamps),1);
    
    parfor s = 1:length(specslope{ii}.timestamps)
        [a b]= min(abs(relTime - specslope{ii}.timestamps(s)));
        [aa bb] = min(abs(abTime(b)-absTime));
        specSlope_t(s) = specslope{ii}.data(s);
        [aa bb] = min(abs(abTime(b)-absTime));
        specSlopeTimes_t(s) = absTime(bb);
        c_slope = 1 + c_slope; 
    end
    specSlopeTimes = [specSlopeTimes;specSlopeTimes_t];
    specSlope = [specSlope;specSlope_t];
    
    slState = SleepState{ii};
    mo = movement{ii};
    rippleAlignedCGM_t = []; ripTime_t =[]; time_t =[]; ripState_t=[];offsets_t=[];
    
    parfor i=1:length(ripples{ii}.peaks)
        [a b]= min(abs(relTime-rips.peaks(i)));
        [aa bb] = min(abs(abTime(b)-absTime));
        
        
        m_t(i) = mo(b);
        zTime_t(i) = zt_rec{ii}(b);
        absolute_t(i) = absTime(bb);
        amplitude_t(i) = rips.data.zScoredAmps(i);
        duration_t(i) = rips.data.duration(i);
        frequency_t(i) = rips.data.peakFrequency(i);
        recording_t(i) = ii;
        thetaPower_z_t(i) = thetaPow_z{ii}(b);
        thetaPower_raw_t(i) = thetaPow_raw{ii}(b);
        
        
        [blah loc] = min(abs(slState.idx.timestamps-rips.peaks(i)));
        state_t(i) = slState.idx.states(loc);
        %% EMG add
        [emg_pk emg_loc]= min(abs(emgFromLFP{ii}.EMGFromLFP.timestamps-rips.peaks(i)));
        emg_idx = emg_loc-75:emg_loc+75;
        emg_idx(emg_idx<1)=[]; emg_idx(emg_idx>length(emgFromLFP{ii}.EMGFromLFP.data))=[]; 
        emg_t(i) = nanmean(emgFromLFP{ii}.EMGFromLFP.data(emg_idx));  % sampled at 0.5 Hz, take avg over 5 min block?? 
        
        
%         [a b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
        
        offset = (abTime(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
        
        idx = find(absTime> abTime(b)-1.15741277113557e-05*60*60*2 & ...
            absTime < abTime(b)+1.15741277113557e-05*60*60*2);
        
        if length(idx)>48 & ~isempty(idx)
            idx = idx(1:48);
        elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
            idx(end+1) = idx(end) + 1;
        end
        if abs(offset) < 60*5 & length(idx) == 48
            rippleAlignedCGM_t(i,:) = isa(idx);
            time_t(i,:) = linspace(-7200,7200,48) - offset;
            ripTime_t(i) = rips.peaks(i);
            if i ~= 1 & i < length(rips.peaks)
                IRI_t(i) = mean(diff(rips.peaks(i-1:i+1)));
            end
            ripState_t(i) = state_t(i);
            offsets_t(i) = offset;
        end

%         m(c) = movement{ii}(b);
%         zTime(c) = zt_rec{ii}(b);
%         absolute(c) = absTime(bb);
%         amplitude(c) = ripples{ii}.data.zScoredAmps(i);
%         duration(c) = ripples{ii}.data.duration(i);
%         frequency(c) = ripples{ii}.data.peakFrequency(i);
%         recording(c) = ii;
%         thetaPower_z(c) = thetaPow_z{ii}(b);
%         thetaPower_raw(c) = thetaPow_raw{ii}(b);
%         
%         
%         [blah loc] = min(abs(SleepState{ii}.idx.timestamps-ripples{ii}.peaks(i)));
%         state(c) = SleepState{ii}.idx.states(loc);
%         %% EMG add
%         [emg_pk emg_loc]= min(abs(emgFromLFP{ii}.EMGFromLFP.timestamps-ripples{ii}.peaks(i)));
%         emg_idx = emg_loc-75:emg_loc+75;
%         emg_idx(emg_idx<1)=[]; emg_idx(emg_idx>length(emgFromLFP{ii}.EMGFromLFP.data))=[]; 
%         emg(c) = nanmean(emgFromLFP{ii}.EMGFromLFP.data(emg_idx));  % sampled at 0.5 Hz, take avg over 5 min block?? 
%         
%         
% %         [a b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
%         
%         offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
%         
%         idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
%             absTime < at{ii}(b)+1.15741277113557e-05*60*60*2);
%         
%         if length(idx)>48 & ~isempty(idx)
%             idx = idx(1:48);
%         elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
%             idx(end+1) = idx(end) + 1;
%         end
%         if abs(offset) < 60*5 & length(idx) == 48
%             rippleAlignedCGM(c,:) = isa(idx);
%             time(c,:) = linspace(-7200,7200,48) - offset;
%             ripTime(c) = ripples{ii}.peaks(i);
%             if i ~= 1 & i < length(ripples{ii}.peaks)
%                 IRI(c) = mean(diff(ripples{ii}.peaks(i-1:i+1)));
%             end
%             ripState(c) = state(c);
%             offsets(c) = offset;
%             c=1+c;
%         end
    end
    
    rippleAlignedCGM = [rippleAlignedCGM; rippleAlignedCGM_t];
%     time = 
    rippleAlignedCGM_t = [];
    ripTime = [ripTime; ripTime_t'];
    ripState = [ripState; ripState_t'];
    offsets = [offsets; offsets_t'];
    
    m = [m, m_t];
    zTime = [zTime, zTime_t];
    absolute = [absolute, absolute_t];
    amplitude = [amplitude,amplitude_t];
    duration = [duration, duration_t];
    frequency = [frequency, frequency_t];
    recording = [recording, recording_t];
    thetaPower_z = [thetaPower_z, thetaPower_z_t];
    thetaPower_raw = [thetaPower_raw, thetaPower_raw_t];
    state = [state, state_t];
    emg = [emg, emg_t];
    clear *_t
    
    ii

end
end
end

count = nan(1,length(absTime));
rec = nan(1,length(absTime));
dur = nan(1,length(absTime));
amps = nan(1,length(absTime));
freq = nan(1,length(absTime));
states = nan(1,length(absTime));
theta_z = nan(1,length(absTime));
theta_raw = nan(1,length(absTime));
spSlope = nan(1,length(absTime));
emgSig = nan(1,length(absTime));

for i=1:length(absTime)
%     ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05*60*5 absTime(i)]));
    ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05 absTime(i)]));
    
    mov(i) = nanmean(m(ind));
    amps(i) = nanmean(amplitude(ind));
    freq(i) = nanmean(frequency(ind));
    dur(i) = nanmean(duration(ind));
%     states(i) = (sum(state(ind)==1) ./ (sum(state(ind)==3))+.0001);  % wake to NREM sleep ratio
    states(i) = (sum(state(ind)==1) - (sum(state(ind)==3)))./(length(ind)+1);  % wake to NREM sleep ratio
    rec(i) = nanmean(recording(ind));
    theta_z(i) = mean(thetaPower_z(ind));
    theta_raw(i) = mean(thetaPower_raw(ind));
    
%     states(i) = nansum(state(ind)==3)./length(ind);
    count(i) = length(ind);
    
    if ~isempty(ind)
        zeitTimes(i) = zTime(ind(end));
    else
        zeitTimes(i) = nan;
    end

    if ~isempty(ind)
        spSlope(i) = nanmean(specSlope(ceil(ind(1)/10):round(ind(end)/10)));
    else
        spSlope(i) = nan;
    end
    emgSig(i) = nanmean(emg(ind));
    
    ind = find(InIntervals(stimTimesAbs,[absTime(i)-1.15741277113557e-05*60*5 absTime(i)]));
    stimRate(i) = length(ind);
    responseMags(i,:) = nanmean(responseMag(ind,:));

end


% idx = find(absTime > 737397.05);   % for CGM1_day1
idx = 1:length(count);
% count(count==0)=nan;clear at zt_rec


% isig_levels(isig_levels<2) = nan;

subplot(3,2,1)
plot(absTime(idx),count(idx),'r')
subplot(3,2,2)
plot(absTime(idx),isig_levels(idx),'k')
subplot(3,2,3)
scatter(count(idx)',isig_levels(idx)','.')
[a b] = corr(count(idx)',isig_levels(idx)','rows','complete');
title(a)
subplot(3,2,4)
for i = -200:200
cc(i+201) = corr(circshift(count(idx),i)',[ 0; diff(isig_levels(idx)')],'rows','complete');
end
plot([-200:200]/5,cc);

subplot(3,2,5)
clear av t
c = diff(count);
for i=-200:200
f = find(c>i&c<i+20);
t = zeros(1,40);
for j=1:length(f)
if f(j) > 20 & f(j) +20<length(isig_levels)
t(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
end
end
av(201+i,:) = nanmean(t,1);clear t;
end
imagesc(av)




figure
for i=1:size(rippleAlignedCGM,1)
time(i,:) = linspace(-7200,7200,48);
time(i,:) = time(i,:) - offsets(i);
end

for i=1:size(rippleAlignedCGM,1)
    r(i,:) = diff(rippleAlignedCGM(i,:));
    rz(i,:) = zscore(rippleAlignedCGM(i,:));
end
timeD = time(:,2:end);
for i= -7200:7200
f = find(timeD(:)>i-120 & timeD(:)<i+120);
avg(i+7201)=nanmean(r(f));
end
plot(avg)



save('/home/david/Dropbox/Documents/pubs/inProgress/glucose/data/cgm8.mat','-v7.3')
% whos
% 
% for i=-60:60
%     pred = [theta',count',amps',dur',freq',states',zt'];
% [beta(i+61,:) dev(i+61) stats(i+61)]=glmfit(pred(idx,:),circshift([0,diff(isig_levels(idx))]',-i),'normal');
% end
% for i=1:121
% p(i,:) = stats(i).p;
% end
% 
% 
% 
% 







