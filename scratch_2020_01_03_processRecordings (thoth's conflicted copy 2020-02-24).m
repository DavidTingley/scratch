
folders = dir('cgm*');
% 
for f = length(folders):-1:1
    cd(folders(f).name)
    
    if exist('amplifier.dat')
       movefile('amplifier.dat',[folders(f).name '.dat']) 
    end
    if exist('amplifier.xml')
       movefile('amplifier.xml',[folders(f).name '.xml']) 
    end
    if exist('amplifier.nrs')
       movefile('amplifier.nrs',[folders(f).name '.nrs']) 
    end
%     
    if ~exist([folders(f).name '.lfp'])
        bz_LFPfromDat(pwd)
        sessionInfo = bz_getSessionInfo(pwd,'noPrompts',true);
        sessionInfo.left = 2;
        sessionInfo.right = 4;
        sessionInfo.ref = 4;
        
        save([folders(f).name '.sessionInfo.mat'],'sessionInfo')
        SleepScoreMaster(pwd)
    end

%     read_Intan_RHD2000_file
%     num_channels = length(boaraddd_adc_channels);  fileinfo = dir('analogin.dat'); num_samples = fileinfo.bytes/(num_channels * 2); fid = fopen('analogin.dat', 'r');
%     v = fread(fid, [num_channels, num_samples], 'uint16'); fclose(fid); v = v * 0.000050354; %
%     v = single(v);
%     [pul] = getPulses_thresh(v,'fs',20000);
%     sessionInfo = bz_getSessionInfo;
%     save([folders(f).name '.stimTimes.mat'],'pul')
% 
%     if ~exist(([folders(f).name '.stimTimes.mat']))
%         fileinfo = dir('digitalin.dat'); 
%         num_samples = fileinfo.bytes/2; % uint16 = 2 bytes 
%         fid = fopen('digitalin.dat', 'r'); 
%         digital_word = fread(fid, num_samples, 'uint16'); 
%         fclose(fid); 
%         ch=0;
%         digital_input_ch = (bitand(digital_word, 2^ch) > 0); % 
%         pul = find(digital_input_ch == 1);
% 
%         c=1;
%         stims = pul(1);
%         for i=2:length(pul)
%             if pul(i)-stims(end) > 201
%             stims(c) = pul(i);
%             c=c+1;
%             end
%         end
%         clear pul;
%         pul{1} = stims;
%         save([folders(f).name '.stimTimes.mat'],'pul')
%     end
    
%     
%     
%     lfp = bz_GetLFP(sessionInfo.channels);
%     
%     [b a]=butter(4,[140/625 180/625],'bandpass');
% 
%     for i=1:length(lfp.channels)
%         filt = FiltFiltM(b,a,single(lfp.data(:,i)));
%         pow = fastrms(filt,15);    
%         mRipple(i) = mean(pow);
%         meRipple(i) = median(pow);
%         mmRippleRatio(i) = mRipple(i)./meRipple(i);
%     endc
% 
%     mmRippleRatio(mRipple<1) = 0;
%     mmRippleRatio(meRipple<1) = 0;validateattributes(Yin,{'double','single'},{'nonempty','real','vector'},...
% z
% 
%     [minVal loc] = max(mmRippleRatio);
%     chan = lfp.channels(loc);
%     plot(mmRippleRatio)
%     
    if ~exist([folders(f).name '.CA1Ripples.events.mat'])
%         SleepScoreMaster(pwd)
        ch = sessionInfo.left; %input('lfp chan: ');
        noi = sessionInfo.ref; %input('noise chan: ');
%         chan = lfp.channels(ch);
%         noise = lfp.channels(noi);
        lfp = bz_GetLFP(ch);
        noise = bz_GetLFP(noi);
        [ripples] = bz_FindRipples(lfp.data,lfp.timestamps,'noise',noise.data,'emgThresh',0,'durations',[20 200]);
        [b a]=butter(4,[140/625 180/625],'bandpass');
        filt = FiltFiltM(b,a,single(lfp.data));
        [ripples.maps,ripples.data,ripples.stats] = bz_RippleStats(double(filt),lfp.timestamps,ripples);
        sessionInfo = bz_getSessionInfo;
        close all
        bz_PlotRippleStats(ripples.maps,ripples.data,ripples.stats)
        close all
        save([sessionInfo.FileName '.CA1Ripples.events.mat'],'ripples')   
    end
    cd ..
end


folders = dir('cgm*');

for f = 1:length(folders)
    cd(folders(f).name)
    if ~exist([folders(f).name '.stimTimes.mat']) & exist('analogin.dat')
        read_Intan_RHD2000_file
        num_channels = length(board_adc_channels);  
        fileinfo = dir('analogin.dat'); 
        num_samples = fileinfo.bytes/(num_channels * 2); 
        fid = fopen('analogin.dat', 'r');
        v = int16(fread(fid, [num_channels, num_samples], 'uint16')); 
        fclose(fid); %v = v * 0.000050354; %
    %     v = single(v);
        [pul] = getPulses_thresh(v,'fs',20000);
        sessionInfo = bz_getSessionInfo;
        save([folders(f).name '.stimTimes.mat'],'pul')
    end
    cd ..
end

[hi lo]=butter(4,[80/625 200/625],'bandpass');

for f = 1:length(folders)
    cd(folders(f).name)
    sessionInfo = bz_getSessionInfo;
    if exist([folders(f).name '.stimTimes.mat']) %& ~exist([sessionInfo.FileName '.responseMagnitudes.mat'])
        
        lfp = bz_GetLFP([sessionInfo.left sessionInfo.right]);
        load([sessionInfo.FileName '.stimTimes.mat'],'pul')
        if length((pul{1})) > 2
            response1 = nan(length(pul{1}),2501);
            response2 = nan(length(pul{1}),2501);
            respMag1 = nan(length(pul{1}),1);
            respMag2 = nan(length(pul{1}),1);
            
        for stim = 1:length(pul{1})
            
            start = ((pul{1}(stim)-1));
            stop = ((pul{1}(stim)+1));
            start = round(start*1250);
            stop = round(stop*1250);
            if stop < length(lfp.timestamps) & start > 0
%             spec = bz_WaveSpec(lfp,'intervals',[start',stop'],'nfreqs',[100],'frange',[1 200],'ncyc',3);
            
%             [specslope] = bz_PowerSpectrumSlope(lfp,3,1/1250,'ints',[start',stop'],'spectype','wavelet','frange',[4 200],'nfreqs',100);

            r1 = makeLength(FiltFiltM(hi,lo,double(lfp.data(start:stop,1))),2501);
            r2 = makeLength(FiltFiltM(hi,lo,double(lfp.data(start:stop,2))),2501);
            
            respMag1(stim) = mean(fastrms(r1(1140:1260),20));
            respMag2(stim) = mean(fastrms(r2(1140:1260),20));
            
            response1(stim,:) = r1;
            response2(stim,:) = r2;
            end 
        end
        response(1,:,:) = response1;
        response(2,:,:) = response2;
        responseMag(:,1) = respMag1;
        responseMag(:,2) = respMag2;
        clear *1 *2
        save([sessionInfo.FileName '.responseMagnitudes.mat'],'res*','-v7.3')
        clear res* spec* lfp pul
        end
    end
    
    cd ..
end














