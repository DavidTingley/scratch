clear all
cd /home/david/Dropbox/Documents/pubs/inProgress/Hormones
% 

animal = 'CGM2'; 

read_CGM_data

cd(['/mnt/packrat/userdirs/david/zpool1/' animal])
% cd /mnt/nyuShare/Buzsakilabspace/david/zpool2/CGM2
d = dir(['*' animal '*']);

for i=1:length(d)
    if d(i).isdir
        cd(d(i).name);
        ripples{i} = bz_LoadEvents(pwd,'ripples');
        if ~isempty(ripples{i})
            [zt_rec{i} rt{i} at{i}] = bz_getZeitgeberTime(d(i).name,20000);
            SleepState{i} = bz_LoadStates(pwd,'SleepState');

            %% add accelerometer data

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
            movement{i} = makeLength(movement{i},f.bytes/6/5000);
            
            fclose(fid);

        %% get stim times
%         if exist('digitalin.dat')
%             fileinfo = dir('digitalin.dat');
%             num_samples = fileinfo.bytes/2;
%             fid = fopen('digitalin.dat','r');
%             dig = fread(fid,num_samples,'uint16=>uint16');
%             digSort = bitand(dig,2^1)>0; clear dig;
%             for chunk = 1:100000:length(digSort)
%                 if chunk + 100000 < length(digSort)
%                     di(chunk:chunk+99999) = int16(diff(digSort(chunk:chunk+100000)));
%                     di(end+1) = digSort(chunk+100000)-digSort(chunk+100001);
%                 else
%                     di(chunk:length(digSort)-1) = int16(diff(digSort(chunk:end)));
%                     di(end+1) = digSort(length(digSort)-1)-digSort(end);
%                 end
%             end
%             clear digSort
%             pos = find(di==1); clear di
%     %         neg = find(d==-1);
%             stimTimes{i} = pos/20000;  % relative timing w/ in recording
%             clear pos
%         else
%            stimTimes{i} = []; 
%         end
        
         %% add power spectrum
%         lfp = bz_GetLFP(24);
%         [specslope{i},spec] = bz_PowerSpectrumSlope(lfp,10,1,'showfig',false);

        else
            specslope{i} = [];
            SleepState{i} = [];
            zt_rec{i} = [];
            rt{i} = [];
            at{i} = [];
        end
        cd(['/mnt/packrat/userdirs/david/zpool1/' animal])
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



c=1;
c_NREM = 1;
c_WAKE = 1;
c_stim = 1;
c_slope = 1;
clear absolute

for ii=1:length(d)
    if  d(ii).isdir
if ~isempty(rt{ii}) & ~isempty(ripples{ii}) & ~isempty(SleepState{ii})
    
%     for s = 1:length(stimTimes{ii})
%         [a b]= min(abs(rt{ii}-stimTimes{ii}(s)));
%         [aa bb] = min(abs(at{ii}(b)-absTime));
%         stimTimesAbs(c_stim) = absTime(bb);
%         c_stim = 1 + c_stim; 
%     end

%     for s = 1:length(specslope{ii}.timestamps)
%         [a b]= min(abs(rt{ii} - specslope{ii}.timestamps(s)));
%         [aa bb] = min(abs(at{ii}(b)-absTime));
%         specSlope(c_slope) = specslope{ii}.data(s);
%         [aa bb] = min(abs(at{ii}(b)-absTime));
%         specSlopeTimes(c_stim) = absTime(bb);
%         c_slope = 1 + c_slope; 
%     end
    
    for i=1:length(ripples{ii}.peaks)
        [a b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
        [aa bb] = min(abs(at{ii}(b)-absTime));
        
        
%         m(c) = movement{ii}(b);
        zTime(c) = zt_rec{ii}(b);
        absolute(c) = absTime(bb);
        amplitude(c) = ripples{ii}.data.peakAmplitude(i);
        duration(c) = ripples{ii}.data.duration(i);
        frequency(c) = ripples{ii}.data.peakFrequency(i);
        
        
        
        
        [blah loc] = min(abs(SleepState{ii}.idx.timestamps-ripples{ii}.peaks(i)));
        state(c) = SleepState{ii}.idx.states(loc);
        
        
        [a b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
        
        offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
        
        idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
            absTime < at{ii}(b)+1.15741277113557e-05*60*60*2);
        
        if length(idx)>48 & ~isempty(idx)
            idx = idx(1:48);
        elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
            idx(end+1) = idx(end) + 1;
        end
        if abs(offset) < 60*5 & length(idx) == 48
            rippleAlignedCGM(c,:) = isig_levels(idx);
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

for i=1:length(absTime)
%     ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05*60*5 absTime(i)]));
    ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05 absTime(i)]));
    
%     mov(i) = nanmean(m(ind));
    amps(i) = nanmean(amplitude(ind));
    freq(i) = nanmean(frequency(ind));
    dur(i) = nanmean(duration(ind));
%     states(i) = (sum(state(ind)==1) ./ (sum(state(ind)==3))+.0001);  % wake to NREM sleep ratio
    states(i) = (sum(state(ind)==1) - (sum(state(ind)==3)))./(length(ind)+1);  % wake to NREM sleep ratio
    
%     states(i) = nansum(state(ind)==3)./length(ind);
    count(i) = length(ind);
    
    if ~isempty(ind)
        zeitTimes(i) = zTime(ind(end));
    else
        zeitTimes(i) = nan;
    end

%     if ~isempty(ind)
%         spSlope(i) = nanmean(specSlope(ind));
%     else
%         spSlope(i) = nan;
%     end
    
%     ind = find(InIntervals(stimTimesAbs,[absTime(i)-1.15741277113557e-05*60*5 absTime(i)]));
%     stimRate(i) = length(ind);


end


idx = find(absTime > 737397.05);   % for CGM1_day1
% idx = 1000:length(count);
count(count==0)=nan;
isig_levels(isig_levels<2) = nan;

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







% 
% isig_levels(3266:end) = nan;
% iisig_levels(1003) = nan;
% isig_levels(1070:1080) = nan;
% isig_levels(1340:1360) = nan;
% isig_levels(1370:1375) = nan;
% isig_levles(1240:1245)=nan;
% isig_levles(1235:1245)=nan;
% isig_levels(1240:1245) = nan;
% isig_levels(2510:2518) = nan;
% isig_levels(3099:3100) = nan;
% isig_levels(3090:3100) = nan;
% 
% whos
% 
% 
% 
% 







