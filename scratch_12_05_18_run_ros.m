clear all
% cd C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose
cd /mnt/nyuShare/Homes/dwt244/glucoseRev/CA1_rec/ros


animal = 'ros'; 
addpath(genpath(pwd))
read_CGM_data_REV
clear a b dat nSteps timeStep ts ISIG_vals dataType dates times % clean up

% cd(['Z:\Homes\dwt244\glucoseRev\ros'])
cd /mnt/nyuShare/Homes/dwt244/glucoseRev/CA1_rec/ros
d = dir(['*' animal '*']);

for i=1:length(d)
    if d(i).isdir
        cd(d(i).name);
        ripples{i} = bz_LoadEvents(pwd,'CA1Ripples');
        sessionInfo = bz_getSessionInfo(pwd,'saveMat',false);
        if ~exist([sessionInfo.FileName '.zeitTimes.mat'])
            [zzt_rec rrt aat] = bz_getZeitgeberTime(d(i).name,20000,20000);
            save([sessionInfo.FileName '.zeitTimes.mat'],'zzt_rec', 'rrt', 'aat');
        else
            load([sessionInfo.FileName '.zeitTimes.mat'],'zzt_rec', 'rrt', 'aat');
        end
        zt_rec{i} = zzt_rec;
        rt{i} = rrt; 
        at{i} = aat;
        cd ..
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
c_spwr = ones(8,1);
clear absolute

for ii= [1:length(d)]
    if  d(ii).isdir
        
if ~isempty(rt{ii}) & ~isempty(ripples{ii})
    
    
%     for j =1:length(ripples{ii}.spwrs)
%         for i=1:length(ripples{ii}.spwrs{j})
%             [a b]= min(abs(rt{ii}-ripples{ii}.spwrs{j}(i)));
%             [aa bb] = min(abs(at{ii}(b)-absTime));
% 
% 
%             spwrTimes{j}(c_spwr(j)) = absTime(bb);
%             c_spwr(j) = c_spwr(j)+1;
%         end
%     end
    for i=1:length(ripples{ii}.peaks)
        [a b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
        [aa bb] = min(abs(at{ii}(b)-absTime));
        
        
        zTime(c) = zt_rec{ii}(b);
        absolute(c) = absTime(bb);
%         amplitude(c) = ripples{ii}.data.zScoredAmps(i);
%         duration(c) = ripples{ii}.data.duration(i);
%         frequency(c) = ripples{ii}.data.peakFrequency(i);
        recording(c) = ii;
     
        
           
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
%             ripState(c) = state(c);
            offsets(c) = offset;
            c=1+c;
        end
    end
    ii
end
end
end


count = nan(1,length(absTime));
rec = nan(1,length(absTime));

for i=1:length(absTime)
    ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
%     ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05 absTime(i)]));
   
    count(i) = length(ind);
    rec(i) = mean(recording(ind));
    
    if ~isempty(ind)
        zeitTimes(i) = zTime(ind(end));
    else
        zeitTimes(i) = nan;
    end
    
%     for j=1:length(ripples{ii}.spwrs)
%        ind = find(InIntervals(spwrTimes{j},[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
%        spwrRates(j,i) = length(ind); 
%     end
end





% save('C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose\data\revision\ros.mat','-v7.3')
save('/home/david/Dropbox/Documents/pubs/inProgress/glucose/data/revision/ros.mat','-v7.3')






