d = dir('*DT12_*');
% ch = [15     7    63    55    47    39    31    23];
ch = 0:63;

% ch = [0    29    12     1];
[b a]=butter(4,[140/625],'high');

for i=1:length(d)
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo;
    try
        check = bz_GetLFP(127,'intervals',[0 1]);
      
        if exist([d(i).name '.CA1Ripples.events.mat'])
            load([d(i).name '.CA1Ripples.events.mat'])
        %     for j=1:4
        %         ch(j) = sessionInfo.spikeGroups.groups{j}(1);
        %     end
            lfp = bz_GetLFP(ch,'intervals',[ripples.peaks-1 ripples.peaks+1]);
            parfor k=1:length(lfp)-1
                for j=1:64
                p(k,j,:) = fastrms(filtfilt(b,a,double(lfp(k).data(:,j))),13);
                end
            end
            s = strsplit(d(i).name,'_');
            s = strsplit(s{2},'u');
            depth(i) = str2num(s{1});
            map(i,:,:) = squeeze(nanmean(p,1));
            for spk = 1:8
                id = sessonInfo.spikeGroups.groups{spk};
                mapAvg(i,spk,:) = nanmean(squeeze(map(i,id+1,:)));
            end
            nRips(i) = length(lfp);
        end
    catch
        disp('nope not this one..')    
    end
    cd ..
end