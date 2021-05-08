d = dir('*DT12_*');
% ch = [15     7    63    55    47    39    31    23];
ch = [11     3    59    51    43    35     27    19];
% ch = 0:63;
% ch = [0    29    12     1];
% [b a]=butter(4,[140/625],'high');
[bb aa]=butter(4,[300/625],'high');
summedWave = zeros(length(d),length(ch),1250,200);

for i=1:length(d)
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo;
%     try
%         check = bz_GetLFP(127,'intervals',[0 1]);
      
        if exist([d(i).name '.CA1Ripples.events.mat']) & exist([sessionInfo.FileName '.lfp'])
            load([d(i).name '.CA1Ripples.events.mat'])
        %     for j=1:4
        %         ch(j) = sessionInfo.spikeGroups.groups{j}(1);
        %     end
%             lfp = bz_GetLFP(ch,'intervals',[ripples.peaks-1 ripples.peaks+1]);
            
%             r = randperm(length(ripples.peaks));
%             if length(ripples.peaks)>50
% %                 lfp = lfp(r(1:500));
%                 ripples.peaks = ripples.peaks(r(1:50));
%             end
            
            for k=1:length(ripples.peaks)-1
                dat = bz_LoadBinary([sessionInfo.FileName '.lfp'],'frequency',1250,'nChannels',sessionInfo.nChannels,'channels',ch+1,'start',ripples.peaks(k)-.5,'duration',1);
            
%                 if length(lfp(k).timestamps) == 2500
                for j=1:length(ch)
%                 p(k,j,:) = fastrmsima(filtfilt(b,a,double(lfp(k).data(:,j))),13);
                    ff = filtfilt(bb,aa,double(dat(:,j)));
                    pdat(k,j,:) = fastrms(ff,20);
                    pdatZ(k,j,:) = zscore(pdat(k,j,11:end-10));
                    [spec] = bz_WaveSpec(dat(:,j),'frange',[1 625],'nfreq',200,'ncyc',3,'space','lin','samplingRate',1250);
%                  wave(k,j,:,:) = abs(spec.data);
                    summedWave(i,j,:,:) = squeeze(summedWave(i,j,:,:)) + abs(spec.data);
                end
%                 p(k,j,:)=nan;
%                 end
            end
            s = strsplit(d(i).name,'_');
            s = strsplit(s{2},'u');
            depth(i) = str2num(s{1});
            avgPow(i,:,:) = squeeze(nanmedian(pdat,1));
            avgPowZ(i,:,:) = squeeze(nanmedian(pdatZ,1));
%             avgWave(i,:,:,:) = squeeze(nanmedian(wave,1));
            clear p pdat wave
%             for spk = 1:8
%                 id = sessionInfo.spikeGroups.groups{spk};
%                 mapAvg(i,spk,:) = nanmean(squeeze(map(i,id+1,:)));
%             end
            nRips(i) = length(ripples.peaks);
            nChan(i) = sessionInfo.nChannels;
            subplot(2,3,1)
            plot(repmat(depth(i),length(ch),1),squeeze(avgPowZ(i,:,625+2-10)),'.k')
            hold on
            plot(repmat(depth(i),length(ch),1),squeeze(avgPowZ(i,:,125)),'.r')
            subplot(2,3,2)
            imagesc(squeeze(mean(avgPowZ)))
            subplot(2,3,3)
            plot(squeeze(mean(mean(avgPowZ))))
            subplot(2,1,2)
            imagesc(zscore(squeeze(mean(mean(summedWave)))',[],2))
            caxis([-3,3])
            pause(.1)
        end
%     catch
%         disp('nope not this one..')    
%     end
    cd ..
end