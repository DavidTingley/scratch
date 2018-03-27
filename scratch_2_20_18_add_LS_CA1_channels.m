d = dir('*201*');
count=1;
for rec = 89:length(d)
    cd(d(rec).name)
    sessionInfo = bz_getSessionInfo
    sessionInfo.FileName
    if ~isfield(sessionInfo,'ls')
        s = input('is there an LS channel: ', 's');
        sessionInfo.ls = str2num(s);
    end
    if ~isfield(sessionInfo,'ca1')
        s = input('is there a CA1 channel: ', 's'); 
        sessionInfo.ca1 = str2num(s);
    end
    if ~isfield(sessionInfo,'refChan')
        s = input('is there a reference channel: ', 's'); 
        sessionInfo.refChan = str2num(s);
    end
    save([sessionInfo.FileName '.sessionInfo.mat'],'sessionInfo')
   cd /home/david/datasets/ripples_LS 
end