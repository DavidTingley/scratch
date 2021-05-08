d = dir('*dt15*');
% ch = [15     7    63    55    47    39    31    23];
ch = 0:31;

% ch = [0    29    12     1];
[b a]=butter(4,[140/625],'high');
summedWave = zeros(length(d),32,2500,200);


for i=1:length(d)
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo(pwd,'noprompts',true);
    try
        check = bz_GetLFP(32,'intervals',[0 1]);
      
        if exist([d(i).name '.CA1Ripples.events.mat'])
            load([d(i).name '.CA1Ripples.events.mat'])
        %     for j=1:4
        %         ch(j) = sessionInfo.spikeGroups.groups{j}(1);
        %     end
            lfp = bz_GetLFP(ch,'intervals',[ripples.peaks-1 ripples.peaks+1]);
            r = randperm(length(lfp));
            if length(r)>1000
                lfp = lfp(r(1:1000));
            end
            for k=1:length(lfp)-1
                for j=1:32
                    p(k,j,:) = fastrms(filtfilt(b,a,double(lfp(k).data(:,j))),13);

                    [spec] = bz_WaveSpec(lfp(k).data(:,j),'frange',[1 400],'nfreqs',200,'ncyc',3,'space','lin','samplingRate',lfp(k).samplingRate);
%                     wave(k,j,:,:) = abs(spec.data);
                    summedWave(i,j,:,:) = squeeze(summedWave(i,j,:,:)) + abs(spec.data);
                end
            end
            s = strsplit(d(i).name,'_');
            s = strsplit(s{2},'u');
            depth(i) = str2num(s{1});
            avgPow(i,:,:) = squeeze(nanmedian(p,1));
%             avgWave(i,:,:,:) = squeeze(nanmedian(wave,1));
            clear p wave
%             for spk = 1:4
%                 id = sessionInfo.spikeGroups.groups{spk};
%                 mapAvg(i,spk,:) = nanmedian(squeeze(map(i,id+1,:)));
%             end
            nRips(i) = length(lfp);
            plot(repmat(depth(i),32,1),squeeze(avgPow(i,:,1251)),'.')
            hold on
            pause(.1)
        end
    catch
        disp('nope not this one..')    
    end
    cd ..
end