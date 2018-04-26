clear all
% cd /home/david/datasets/ripples_LS
recordingList  = dir('*201*');
animals = {'DT1','DT2','DT5','DT7','DT8','DT9'};

for i=1:length(recordingList)
    cd(recordingList(i).name);
    spikes = bz_GetSpikes;
    if ~isempty(spikes)
    for spk = 1:length(spikes.times)
        if strcmp(spikes.region{spk},'ls')
            reg(spikes.spindices(:,2)==spk) = 1;
        end
        if strcmp(spikes.region{spk},'hpc') | strcmp(spikes.region{spk},'ca1') |strcmp(spikes.region{spk},'ca3')
            reg(spikes.spindices(:,2)==spk) = 2;
        end
    end
    if length(unique(reg))>1 & sum(unique(reg))==3 & length(reg)==length(spikes.spindices)
        reg(reg==0) = 3;
        [ccg t] = CCG(spikes.spindices(:,1),reg,'binSize',.001,'duration',1);
        cc(:,i) = ccg(:,1,2);
        nCells(i) = spk;
    end
    clear reg
    end
    cd ..    
end


