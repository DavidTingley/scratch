addpath(genpath('/home/david/Dropbox/code'))
cd /home/david/datasets/ripples_LS/
% cd E:\datasets\ripples_LS
clear all
ripReward = nan(98,201);

d = dir('*201*');

for rec = 1:length(d)
    cd(d(rec).name)
    if exist(([d(rec).name '.behavior.mat']))
        sessionInfo = bz_getSessionInfo;
        sessionInfo.FileName;
        ca1_ripples = bz_LoadEvents(pwd,'CA1Ripples_4SD');
        
        ls_ripples = bz_LoadEvents(pwd,'LSRipples');
        
        load([d(rec).name '.behavior.mat']);
        if ~isempty(ca1_ripples) & ~isempty(behavior)
            [times groups] = spikes2sorted({behavior.events.trialIntervals(:,2),ca1_ripples.peaks});
            [ccg t] = CCG(times,groups,'binSize',.01,'duration',2,'norm','rate');
            ripReward(rec,:) = ccg(:,1,2);
            imagesc(ripReward)
            pause(.1)
            
            
            for trial = 1:length(behavior.events.trialIntervals(:,2))
                nRips{rec}(trial) = length(Restrict(ca1_ripcd ples.peaks,[behavior.events.trialIntervals(trial,2) behavior.events.trialIntervals(trial,2)+.75]));
            end 
        end
        
        
    end
   cd /home/david/datasets/ripples_LS 

end

