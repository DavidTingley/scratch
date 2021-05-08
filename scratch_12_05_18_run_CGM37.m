clear all
% for threshold = [1000 2000 5000 7500 10000 15000 20000 25000 30000 40000 50000 75000 100000]
%     if ~exist(['/home/david/Dropbox/Documents/pubs/inProgress/glucose/data/revision/dvSweep/' num2str(threshold,'%06.f') '_Vanessa.mat'])

% cd C:\Users\SB13FLLT001\Dropbox\Documents\pubs\inProgress\glucose
cd /mnt/nyuShare/Homes/dwt244/glucoseRev/dvHPC_rec/CGM37

animal = 'CGM37'; 
addpath(genpath(pwd))
read_CGM_data_REV
clear a b dat nSteps timeStep ts ISIG_vals dataType dates times % clean up

% cd(['Z:\Homes\dwt244\glucoseRev\bruce'])
cd  /mnt/nyuShare/Homes/dwt244/glucoseRev/dvHPC_rec/CGM37
d = dir(['*' animal '*']);

for i=1:length(d)
    if d(i).isdir
        cd(d(i).name);
%         dRipples{i} = bz_LoadEvents(pwd,[num2str(threshold,'%06.f') '.DorsalRipples']);
%         vRipples{i} = bz_LoadEvents(pwd,[num2str(threshold,'%06.f') '.VentralRipples']);
        dRipples{i} = bz_LoadEvents(pwd,['DorsalRipples']);
        vRipples{i} = bz_LoadEvents(pwd,['VentralRipples']);
        gRipples{i} = bz_LoadEvents(pwd,['GlobalRipples']);
        vsharpWaves{i} = bz_LoadEvents(pwd,'sharpWaves');
        sessionInfo = bz_getSessionInfo;
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


remove = [1:72 632:637];
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
vc = 1;
vg = 1;
vsc = 1;
c_NREM = 1;
c_WAKE = 1;
c_stim = 1;
c_slope = 1;
c_mua = 1; 
clear absolute vabsolute vabsolute_sw

for ii= [1:length(d)]
    if  d(ii).isdir
        % ventral sharpwaves
       if ~isempty(rt{ii}) & ~isempty(vsharpWaves{ii})
        for i=1:length(vsharpWaves{ii}.timestamps)
            [a b]= min(abs(rt{ii}-vsharpWaves{ii}.timestamps(i)));
            [aa bb] = min(abs(at{ii}(b)-absTime));

            vabsolute_sw(vsc) = absTime(bb);
            
            offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
        
            idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
                absTime < at{ii}(b)+1.15741277113557e-05*60*60*2);

            if length(idx)>48 & ~isempty(idx)
                idx = idx(1:48);
            elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
                idx(end+1) = idx(end) + 1;
            end
            if abs(offset) < 60*5 & length(idx) == 48
                vsc=1+vsc;
            end
        
        end
       end
    
    if ~isempty(rt{ii}) & ~isempty(gRipples{ii})
        for i=1:length(gRipples{ii}.peaks)
            [a b]= min(abs(rt{ii}-gRipples{ii}.peaks(i)));
            [aa bb] = min(abs(at{ii}(b)-absTime));

            gabsolute(vg) = absTime(bb);
            
            offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
        
            idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
                absTime < at{ii}(b)+1.15741277113557e-05*60*60*2);

            if length(idx)>48 & ~isempty(idx)
                idx = idx(1:48);
            elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
                idx(end+1) = idx(end) + 1;
            end
            if abs(offset) < 60*5 & length(idx) == 48
                vg=1+vg;
            end
        
        end
    end
    
    if ~isempty(rt{ii}) & ~isempty(vRipples{ii})
        for i=1:length(vRipples{ii}.peaks)
            [a b]= min(abs(rt{ii}-vRipples{ii}.peaks(i)));
            [aa bb] = min(abs(at{ii}(b)-absTime));

            vabsolute(vc) = absTime(bb);
            
            offset = (at{ii}(b) - absTime(bb)) / 1.15741277113557e-05;  % offset in seconds
        
            idx = find(absTime> at{ii}(b)-1.15741277113557e-05*60*60*2 & ...
                absTime < at{ii}(b)+1.15741277113557e-05*60*60*2);

            if length(idx)>48 & ~isempty(idx)
                idx = idx(1:48);
            elseif length(idx)< 48  & ~isempty(idx) & idx(end) ~= length(absTime)
                idx(end+1) = idx(end) + 1;
            end
            if abs(offset) < 60*5 & length(idx) == 48
                vc=1+vc;
            end
        
        end
    end
    
if ~isempty(rt{ii}) & ~isempty(dRipples{ii})
    
    
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

    for i=1:length(dRipples{ii}.peaks)
        [a b]= min(abs(rt{ii}-dRipples{ii}.peaks(i)));
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
            if i ~= 1 & i < length(dRipples{ii}.peaks)
                IRI(c) = mean(diff(dRipples{ii}.peaks(i-1:i+1)));
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
vcount = nan(1,length(absTime));
vcount_sw = nan(1,length(absTime));
rec = nan(1,length(absTime));

for i=1:length(absTime)
    if exist('absolute')
        ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
    %     ind = find(InIntervals(absolute,[absTime(i)-1.15741277113557e-05 absTime(i)]));

        count(i) = length(ind);
        rec(i) = mean(recording(ind));

        if ~isempty(ind)
            zeitTimes(i) = zTime(ind(end));
        else
            zeitTimes(i) = nan;
        end
    end
    
    
    if exist('gabsolute')
       ind = find(InIntervals(vabsolute,[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
       gcount(i) = length(ind); 
    end
    
    if exist('vabsolute')
       ind = find(InIntervals(vabsolute,[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
       vcount(i) = length(ind); 
    end
   
    if exist('vabsolute_sw')
       ind = find(InIntervals(vabsolute_sw,[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
       vcount_sw(i) = length(ind); 
    end
%     for j=1:length(ripples{ii}.spwrs)
%        ind = find(InIntervals(spwrTimes{j},[absTime(i)-1.15741277113557e-05*60*2.5 absTime(i)+1.15741277113557e-05*60*2.5]));
%        spwrRates(j,i) = length(ind); 
%     end
end

count(isnan(rec))=nan;
gcount(isnan(rec))=nan;
vcount(isnan(rec))=nan;
vcount_sw(isnan(rec))=nan;
idx = 45:637;

clear dRipples vRipples
% save(['/home/david/Dropbox/Documents/pubs/inProgress/glucose/data/revision/dvSweep/' num2str(threshold,'%06.f') '_CGM37.mat'],'-v7.3')
save(['/home/david/Dropbox/Documents/pubs/inProgress/glucose/data/revision/'  'CGM37.mat'],'-v7.3')
%     end
% end



