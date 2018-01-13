d = dir('*201*');

for rec = 1:length(d)
    cd(d(i).name)
    sessionInfo = bz_getSessionInfo;
    if ~isempty(sessionInfo.ca1)
        ca1 = bz_getSessionInfo.ca1;
        ls = bz_getSessionInfo.ls;
        
        [ripples]=bz_FindRipples(ca1.data,ca1.timestamps),
        
        
        
    end
end