d = dir('*201*');

for ii=1:length(d)
    cd(d(ii).name)
    load([d(ii).name '.behavior.mat'])
    if strcmp(behavior.units, 'pixels')
        
%%
    b = dir('*pos');
    pos = importdata(b.name);
%     for ii=1:5
%         p(:,ii) = makeLength(pos(:,ii),length(pos)*4);
%     end
    [interpolatedBehav,trials] = bz_GetBehavFromPos(behavior, p);
    save([d(i).name '.interpolatedBehav.behavior.mat'],'interpolatedBehav')

%%
%     load([d(i).name '.interpolatedBehav.behavior.mat'],'interpolatedBehav')
%     spikes = bz_GetSpikes;
%     sessionInfo = bz_getSessionInfo;
%     lfp = bz_GetLFP(sessionInfo.thetaChans(end));
%     [firingMaps] = bz_firingMap1D(spikes,interpolatedBehav,lfp,5);
%     save([d(i).name '.firingMaps.cellinfo.mat'],'firingMaps')

%%
% scratch_07_11_17_spatialBinnedGLM

%%
% disp(['starting ' d(ii).name])
% wrapInfoAnalysis
% clearvars -except d ii
    end
    cd /home/david/datasets/lsDataset
    clear p pos b behavior 
end