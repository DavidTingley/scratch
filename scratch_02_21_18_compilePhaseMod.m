d = dir('*201*');
count=1;
vals = [];
for rec = 1:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo;
    if ~isfield(sessionInfo,'animal')
        sessionInfo
        pwd
        s = input('which animal is this? ','s');
        sessionInfo.animal = s;
        save([sessionInfo.FileName '.sessionInfo.mat'],'sessionInfo')
    end
if isfield(sessionInfo,'animal')
    animal = str2num(sessionInfo.animal(3:end));
%     pax = subplot(4,3,animal,polaraxes('rlim',[0 1]));
     subplot(4,3,animal)
    spikes = bz_GetSpikes('noprompts',true,'savemat',false);
    if exist([sessionInfo.FileName '.PhaseLocking_LS_rippleband.mat'])
        load([sessionInfo.FileName '.PhaseLocking_LS_rippleband.mat'],'PhaseLockingData');
        for sp = 1:length(PhaseLockingData)
            uid = PhaseLockingData{sp}.UID;
            ma = max(spikes.rawWaveform{find(uid==spikes.UID)});
            mi = min(spikes.rawWaveform{find(uid==spikes.UID)});
           if abs(ma) < abs(mi) & strcmp(PhaseLockingData{sp}.region{1},'ls')% |strcmp(PhaseLockingData{sp}.region{1},'ca3')|strcmp(PhaseLockingData{sp}.region{1},'ca1') 
           if length(PhaseLockingData{sp}.spkphases{1})  > 50
%                subplot(3,3,spikes.shankID(sp))
               polarplot(PhaseLockingData{sp}.phasestats.m,PhaseLockingData{sp}.phasestats.r,'.k')
               vals = [vals;animal,PhaseLockingData{sp}.phasestats.p,PhaseLockingData{sp}.phasestats.m,PhaseLockingData{sp}.phasestats.r];
               rlim([0 1])
               hold on
           end
           end
        end
%         title([sessionInfo.FileName])
%         pause
%         clf
    end
end
   cd /home/david/datasets/ripples_LS 
   pause(.01)
end


