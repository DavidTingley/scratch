clear all
% cd /home/david/Dropbox/Documents/pubs/inProgress/glucose
cd C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose

animal = 'dt15'; 
addpath(genpath(pwd))
read_CGM_data
clear a b dat nSteps timeStep ts ISIG_vals dataType dates times % clean up

% cd(['/mnt/packrat/userdirs/david/zpool4/hypthal/dt15'])
cd(['D:\datasets\glucose\dt15\'])
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
            if 1 %~exist([d(i).name '.movement.mat']) 
                num_channels = 6;
%                 filePath = ['Z:\dwt244\zpool4\hypthal\dt15\' d(i).name '\'];
                filePath = '';[d(i).name '\'];
                fid = fopen([filePath 'auxiliary.dat'], 'r');
                f = dir([filePath 'auxiliary.dat']);
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
            elseif exist([d(i).name '.movement.mat'])
                load([d(i).name '.movement.mat'],'movement_s')
                movement{i} = movement_s; clear movement_s;
            end
            
%             %% get stim times
%             if ~exist([d(i).name '.stimTimes.mat'])
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
%             else
%                 load([d(i).name '.stimTimes.mat'],'stimTimes_s')
%                 stimTimes{i} = stimTimes_s; clear stimTimes_s
%             end
        
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
            if ~exist([d(i).name '.thetaResid.mat'])
                lfp.data = ripples{i}.detectorinfo.detectionparms.lfp;
                lfp.timestamps = [0.00001:length(lfp.data)]'./1250;
                lfp.samplingRate = 1250;
                [specslope_singleRec,spec] = bz_PowerSpectrumSlope(lfp,20,10,'showfig',false,'nfreqs',100);
                id = find(specslope_singleRec.freqs>5 & specslope_singleRec.freqs<13);
                r = max(specslope_singleRec.resid(:,id)');
                theta_resid = makeLength(fastrms(r,30),ceil(length(lfp.data)./1250));
                save([d(i).name '.thetaResid.mat'],'theta_resid')
                
                thetaPow_resid{i} = theta_resid;  clear theta_resid
            else
                load([d(i).name '.thetaResid.mat'],'theta_resid')
                thetaPow_resid{i} = theta_resid;  clear theta_resid
%                 thetaPow_z{i} = theta_resid;
            end
            if exist([d(i).name '.thetaPower.mat'])
                [bb aa] = butter(3,[6/625 11/625],'bandpass');
                temp = downsample(filtfilt(bb,aa,double(ripples{i}.detectorinfo.detectionparms.lfp)),1250);
%                 temp(1:10)=nan; temp(end-10:end)=nan; % kill those edges
                thetaPow_raw_s = fastrms(temp,60*5);
                thetaPow_z_s = nanZscore(thetaPow_raw_s);
                save([d(i).name '.thetaPower.mat'],'thetaPow*')
%                 thetaPow_raw{i} = thetaPow_raw_s; clear thetaPow_raw_s
                thetaPow_z{i} = thetaPow_z_s; clear thetaPow_z_s
            else
                load([d(i).name '.thetaPower.mat'],'thetaPow_z_s')
%                 thetaPow_raw{i} = thetaPow_raw_s; clear thetaPow_raw_s
                thetaPow_z{i} = thetaPow_z_s; clear thetaPow_z_s
            end
            
            
        else
            emgFromLFP{i} = [];
            thetaPow_z{i} = [];
            thetaPow_resid{i} = [];
            specslope{i} = [];
            SleepState{i} = [];
            zt_rec{i} = [];
            rt{i} = [];
            at{i} = [];
        end
%         cd(['/mnt/packrat/userdirs/david/zpool4/hypthal/' animal])
        cd(['D:\datasets\glucose\dt15\'])
        
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



c=1;
c_NREM = 1;
c_WAKE = 1;
c_stim = 1;
c_slope = 1;
c_mua = 1; 
clear absolute

for ii= [1:length(d)]
    if  d(ii).isdir
if ~isempty(rt{ii}) & ~isempty(ripples{ii})
    
%     for s = 1:length(stimTimes{ii})
%         [a b]= min(abs(rt{ii}-stimTimes{ii}(s)));
%         [aa bb] = min(abs(at{ii}(b)-absTime));
%         stimTimesAbs(c_stim) = absTime(bb);
%         c_stim = 1 + c_stim; 
%     end

%     for s = 1:length(pow{ii}.data)
%         [a b] = min(abs(rt{ii} - pow{ii}.timestamps(s)));
% %         [aa bb] = min(abs(at{ii}(b)-absTime));
%         muaPower(c_mua) = pow{ii}.data(b);
%         [aa bb] = min(abs(at{ii}(b)-absTime));
%         muaTimes(c_mua) = absTime(bb);
%         c_mua = 1 + c_mua; 
%     end

    for s = 1:length(specslope{ii}.timestamps)
        [a b]= min(abs(rt{ii} - specslope{ii}.timestamps(s)));
        [aa bb] = min(abs(at{ii}(b)-absTime));
        specSlope(c_slope) = specslope{ii}.data(s);
        specSlopeTimes(c_slope) = absTime(bb);
        c_slope = 1 + c_slope; 
    end
    
    for i=1:length(ripples{ii}.peaks)
        [a b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
        [aa bb] = min(abs(at{ii}(b)-absTime));
        
        
        m(c) = movement{ii}(b);
        zTime(c) = zt_rec{ii}(b);
        absolute(c) = absTime(bb);
        amplitude(c) = ripples{ii}.data.zScoredAmps(i);
        duration(c) = ripples{ii}.data.duration(i);
        frequency(c) = ripples{ii}.data.peakFrequency(i);
        recording(c) = ii;
        thetaPower_z(c) = thetaPow_z{ii}(b);
        thetaPower_resid(c) = thetaPow_resid{ii}(b);
        
        
        
        [blah loc] = min(abs(SleepState{ii}.idx.timestamps-ripples{ii}.peaks(i)));
        state(c) = SleepState{ii}.idx.states(loc);
        
        %% EMG add
        [emg_pk emg_loc]= min(abs(emgFromLFP{ii}.EMGFromLFP.timestamps-ripples{ii}.peaks(i)));
        emg_idx = emg_loc-75:emg_loc+75;
        emg_idx(emg_idx<1)=[]; emg_idx(emg_idx>length(emgFromLFP{ii}.EMGFromLFP.data))=[]; 
        emg(c) = nanmean(emgFromLFP{ii}.EMGFromLFP.data(emg_idx));  % sampled at 0.5 Hz, take avg over 5 min block?? 
        
%         [a b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
        
        offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
        
        idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
            absTime < at{ii}(b)+1.15741277113557e-05*60*60*2);
        
        if length(idx)>48 & ~isempty(idx)
            idx = idx(1:48);
        elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
            idx(end+1) = idx(end) + 1;
        end
        if abs(offset) < 60*5 & length(idx) == 48
            rippleAlignedCGM(c,:) = isa(idx);
            time(c,:) = linspace(-7200,7200,48) - offset;
            ripTime(c) = absTime(bb);
            if i ~= 1 & i < length(ripples{ii}.peaks)
                IRI(c) = mean(diff(ripples{ii}.peaks(i-1:i+1)));
            end
            ripState(c) = state(c);
            offsets(c) = offset;
            c=1+c;
        end
    end
    ii

%     rips_NREM = Restrict(ripples{ii}.peaks,double(SleepState{ii}.ints.NREMstate));
%     for i=1:length(rips_NREM) 
%         [a b]= min(abs(rt{ii}-rips_NREM(i)));
%         [aa bb] = min(abs(at{ii}(b)-absTime));
%         offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
%         idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
%             absTime< at{ii}(b)+1.15741277113557e-05*60*60*2);
%         if length(idx)>48 & ~isempty(idx)
%             idx = idx(1:48);
%         elseif length(idx)< 48 & ~isempty(idx)
%             idx(end+1) = idx(end) + 1;
%         end1
%         if abs(offset) < 60*5 & length(idx) == 48
%             rippleAlignedCGM_NREM(c_NREM,:) = circshift(makeLength(isig_levels(idx),length(idx)*5),round(offset));
%             c_NREM=1+c_NREM;
%         end
%     end
% 
%     rips_WAKE = Restrict(ripples{ii}.peaks,double(SleepState{ii}.ints.WAKEstate));
%     for i=1:length(rips_WAKE)
%         [a b]= min(abs(rt{ii}-rips_WAKE(i)));
%         [aa bb] = min(abs(at{ii}(b)-absTime));
%         offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
%         idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
%             absTime< at{ii}(b)+1.15741277113557e-05*60*60*2);
%         if length(idx)>48 & ~isempty(idx)
%             idx = idx(1:48);
%         elseif length(idx)< 48 & ~isempty(idx)
%             idx(end+1) = idx(end) + 1;
%         end
%         if abs(offset) < 60*5 & length(idx) == 48
%             rippleAlignedCGM_WAKE(c_WAKE,:) = circshift(makeLength(isig_levels(idx),length(idx)*5),round(offset));
%             c_WAKE=1+c_WAKE;
%         end
%     end
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
theta_resid = nan(1,length(absTime));
spSlope = nan(1,length(absTime));
emgSig = nan(1,length(absTime));

for i=1:length(absTime)
    ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
%     ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05 absTime(i)]));
%     
    mov(i) = nanmean(m(ind));
    amps(i) = nanmean(amplitude(ind));
    freq(i) = nanmean(frequency(ind));
    dur(i) = nanmean(duration(ind));
%     states(i) = (sum(state(ind)==1) ./ (sum(state(ind)==3))+.0001);  % wake to NREM sleep ratio
    states(i) = (sum(state(ind)==1) - (sum(state(ind)==3)))./(length(ind)+1);  % wake to NREM sleep ratio
    
%     states(i) = nansum(state(ind)==3)./length(ind);
    count(i) = length(ind);
    rec(i) = mean(recording(ind));
    theta_z(i) = mean(thetaPower_z(ind));
    theta_resid(i) = mean(thetaPower_resid(ind));
    emgSig(i) = nanmean(emg(ind));
    if ~isempty(ind)
        zeitTimes(i) = zTime(ind(end));
    else
        zeitTimes(i) = nan;
    end
    
    ind = find(InIntervals(specSlopeTimes,[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
    if ~isempty(ind)
        spSlope(i) =  nanmean(specSlope(ind));
    else
        spSlope(i) = nan;
    end
    
%     ind = find(InIntervals(muaTimes,[absTime(i)-1.15741277113557e-05*60*5 absTime(i)]));
%     stimRate(i) = length(ind);
%     MUA_power(i) = mean(muaPower(ind));

end


% idx = find(absTime > 737397.05);   % for CGM1_day1
idx = [1:length(count)]
% idx = 1000:length(count);
count(count==0)=nan;
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
imagesc(-20:20,-200:200,av)




%%%%%%%%%%%%

for i=1:size(rippleAlignedCGM,1)
time(i,:) = linspace(-7200,7200,48);
time(i,:) = time(i,:) - offsets(i);
end
timeD = time(:,2:end);
clear r

for i=1:size(rippleAlignedCGM,1)
r(i,:) = diff(rippleAlignedCGM(i,:));
end
timeD = time(:,2:end);
for i= -7200:7200
f = find(timeD(:)>i-120 & timeD(:)<i+120);
avg(i+7201)=nanmean(r(f));
end
plot(avg)



save('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data\_dt15_latest.mat','-v7.3')
% whos
% 
% 
% 
% 







