clear all
% cd /home/david/Dropbox/Documents/pubs/inProgress/glucose
% cd C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose

% parpool(4)
animal = 'CGM36'; 
addpath(genpath(pwd));
read_CGM_data_REV
clear a b dat nSteps timeStep ts ISIG_vals datatype dates times % clean up


% read_stimLog_data
% stimTimes_logFile = stimTimes_logFile(absTime(end)>stimTimes_logFile);

cd(['/mnt/nyuShare/Homes/dwt244/glucoseRev/CA1_DLX/' animal])
% cd(['D:\datasets\glucose\' animal])
% cd /mnt/nyuShare/Buzsakilabspace/david/zpool2/CGM2
d = dir(['*' animal '*']);

for i=1:length(d)
    if d(i).isdir
        cd(d(i).name);

            %% load some stuff
            if exist([d(i).name '.lfp'])
            [zt_rec{i} rt{i} at{i}] = bz_getZeitgeberTime(d(i).name,20000,20000);
            else
            zt_rec{i} = [];
            rt{i} = [];
            at{i} = [];
            end
            
            %% add accelerometer data
          
            
            %% get stim times0;
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
%                 load([d(i).name '.responseMagnitudes.mat'])
                load([d(i).name '.stimTimes.mat'],'stim*')

                stimTimes1{i} = stim1;
                stimTimes2{i} = stim2; clear stim2 stim1
            end
        
             %% add power spectrum
%             if ~exist([d(i).name '.specSlope.mat'])
%                 lfp.data = ripples{i}.detectorinfo.detectionparms.lfp;
%                 lfp.timestamps = [0.00001:length(lfp.data)]'./1250;
%                 lfp.samplingRate = 1250;
%                 [specslope_singleRec,spec] = bz_PowerSpectrumSlope(lfp,20,10,'showfig',false,'nfreqs',10);
%                 save([d(i).name '.specSlope.mat'],'specslope_singleRec','spec')
%                 specslope{i} = specslope_singleRec; clear specslope_singleRec;
%             else
%                 load([d(i).name '.specSlope.mat'])
%                 specslope{i} = specslope_singleRec; clear specslope_singleRec;
%             end

            %% theta power
%             if ~exist([d(i).name '.thetaPower.mat'])
%                 [bb aa] = butter(3,[6/625 11/625],'bandpass');
%                 temp = decimate(filtfilt(bb,aa,double(ripples{i}.detectorinfo.detectionparms.lfp)),1250);
%                 temp(1:10)=nan; temp(end-10:end)=nan; % kill those edges
%                 thetaPow_raw_s = fastrms(temp,60*5);
%                 thetaPow_z_s = nanZscore(thetaPow_raw_s);
%                 save([d(i).name '.thetaPower.mat'],'thetaPow_raw_s','thetaPow_z_s')
%                 thetaPow_raw{i} = thetaPow_raw_s; clear thetaPow_raw_s
%                 thetaPow_z{i} = thetaPow_z_s; clear thetaPow_z_s
%             else
%                 load([d(i).name '.thetaPower.mat'],'thetaPow_raw_s','thetaPow_z_s')
%                 thetaPow_raw{i} = thetaPow_raw_s; clear thetaPow_raw_s
%                 thetaPow_z{i} = thetaPow_z_s; clear thetaPow_z_s
%             end
                
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
        cd(['/mnt/nyuShare/Homes/dwt244/glucoseRev/CA1_DLX/' animal])
%         cd(['D:\datasets\glucose\' animal])
%         cd /mnt/nyuShare/Buzsakilabspace/david/zpool2/CGM1
       
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
stimTime =[];
stimTimesAbs = [];
stimAlignedCGM2 = [];
stimTime2 = [];
stimTimesAbs2 = [];

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
if ~isempty(rt{ii}) & ~isempty(stimTimes1{ii}) % & ~isempty(SleepState{ii})
    relTime = rt{ii}; 
    abTime = at{ii}; 
    
    sti = stimTimes1{ii};
    if ~isempty(stimTimes1{ii})
    stimAlignedCGM_t = nan(length(stimTimes1{ii}),48);
%     responseMag_t = nan(length(stimTimes{ii}),1251);
    stimTime_t = nan(length(stimTimes1{ii}),1);
    stimTimesAbs_t = nan(length(stimTimes1{ii}),1);
    
    for s = 1:length(sti)
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
            stimAlignedCGM_t(s,:) = isa(idx);
            stimTime_t(s) = abTime(b);
        else
            stimAlignedCGM_t(s,:) = nan(48,1);
            stimTime_t(s) = abTime(b);
        end
    end
    stimAlignedCGM = [stimAlignedCGM;stimAlignedCGM_t];
%     responseMagnitude = [responseMagnitude;responseMag_t];
    stimTime = [stimTime;stimTime_t];
    stimTimesAbs = [stimTimesAbs;stimTimesAbs_t];
    end
    
    %% and two
    sti = stimTimes2{ii};
    if ~isempty(stimTimes2{ii})
    stimAlignedCGM_t = nan(length(stimTimes2{ii}),48);
%     responseMag_t = nan(length(stimTimes{ii}),1251);
    stimTime_t = nan(length(stimTimes2{ii}),1);
    stimTimesAbs_t = nan(length(stimTimes2{ii}),1);
    
    for s = 1:length(sti)
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
            stimAlignedCGM_t(s,:) = isa(idx);
            stimTime_t(s) = abTime(b);
        else
            stimAlignedCGM_t(s,:) = nan(48,1);
            stimTime_t(s) = abTime(b);
        end
    end
    stimAlignedCGM2 = [stimAlignedCGM2;stimAlignedCGM_t];
%     responseMagnitude = [responseMagnitude;responseMag_t];
    stimTime2 = [stimTime2;stimTime_t];
    stimTimesAbs2 = [stimTimesAbs2;stimTimesAbs_t];
    end
    
%%
  
    
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
%     ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05 absTime(i)]));
    
%     mov(i) = nanmean(m(ind));
%     amps(i) = nanmean(amplitude(ind));
%     freq(i) = nanmean(frequency(ind));
%     dur(i) = nanmean(duration(ind));
% %     states(i) = (sum(state(ind)==1) ./ (sum(state(ind)==3))+.0001);  % wake to NREM sleep ratio
%     states(i) = (sum(state(ind)==1) - (sum(state(ind)==3)))./(length(ind)+1);  % wake to NREM sleep ratio
%     rec(i) = nanmean(recording(ind));
%     theta_z(i) = mean(thetaPower_z(ind));
%     theta_raw(i) = mean(thetaPower_raw(ind));
%     
% %     states(i) = nansum(state(ind)==3)./length(ind);
%     count(i) = length(ind);
    
%     if ~isempty(ind)
%         zeitTimes(i) = zTime(ind(end));
%     else
%         zeitTimes(i) = nan;
%     end

%     if ~isempty(ind)
%         spSlope(i) = nanmean(specSlope(ceil(ind(1)/10):round(ind(end)/10)));
%     else
%         spSlope(i) = nan;
%     end
%     emgSig(i) = nanmean(emg(ind));
    
    ind = find(InIntervals(stimTimesAbs,[absTime(i)-1.15741277113557e-05*60*5 absTime(i)]));
    stimRate(i) = length(ind);
    
    ind = find(InIntervals(stimTimesAbs2,[absTime(i)-1.15741277113557e-05*60*5 absTime(i)]));
    stimRate2(i) = length(ind);

end




save('/media/data/Dropbox/Documents/pubs/inProgress/glucose/data/revision/CGM36.mat','-v7.3')






