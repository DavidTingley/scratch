d = dir('*dt15*');
% ch = [15     7    63    55    47    39    31    23];
ch = 0:31;

% ch = [0    29    12     1];
[b a]=butter(4,[140/625],'high');

for i=1:length(d)
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo;
    try
        check = bz_GetLFP(32,'intervals',[0 1]);
      
        if exist([d(i).name '.CA1Ripples.events.mat'])
            load([d(i).name '.CA1Ripples.events.mat'])
        %     for j=1:4
        %         ch(j) = sessionInfo.spikeGroups.groups{j}(1);
        %     end
            lfp = bz_GetLFP(ch,'intervals',[ripples.peaks-1 ripples.peaks+1]);
            for k=1:length(lfp)-1
                for j=1:32
                p(k,j,:) = fastrms(filtfilt(b,a,double(lfp(k).data(:,j))),13);
                end
            end
            s = strsplit(d(i).name,'_');
            s = strsplit(s{2},'u');
            depth(i) = str2num(s{1});
            map(i,:,:) = squeeze(nanmedian(p,1));
            for spk = 1:4
                id = sessionInfo.spikeGroups.groups{spk};
                mapAvg(i,spk,:) = nanmedian(squeeze(map(i,id+1,:)));
            end
            nRips(i) = length(lfp);
        end
    catch
        disp('nope not this one..')    
    end
    cd ..
end