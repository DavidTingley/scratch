d = dir('*201*');
c=1;
region = [];
proportion =[];
state = [];
state_ls = [];
ls_prop = [];

for i=1:length(d)
    cd(d(i).name)
    
    rips = bz_LoadEvents(pwd,'CA1Ripples');
%        rips = bz_LoadEvents(pwd,'popBursts');

    hpc = bz_GetSpikes('noprompts',true,'region','hpc');
    
    ls = bz_GetSpikes('noprompts',true,'region','ls');
    
    if isempty(hpc)
        hpc =bz_GetSpikes('noprompts',true,'region','ca3');
        
    end
    SleepState = bz_LoadStates(pwd,'SleepState');
    
    if ~isempty(hpc) & ...
            ~isempty(rips) & ...
            isfield(SleepState.ints,'NREMstate')
        nCells = length(hpc.times);
        [sleep] = double(InIntervals(rips.peaks,SleepState.ints.NREMstate));
        for r = 1:length(rips.timestamps)
            spks = (InIntervals(hpc.spindices(:,1),[rips.timestamps(r,1)-.025 rips.timestamps(r,2)+.025]));
            nSpks(r) = length(unique(hpc.spindices(spks,2)));
            reg(r) = sum(double(hpc.region{1}));
        end
        state = [state;sleep];
        region = [region, reg];
        proportion = [proportion, nSpks./nCells]; clear nSpks reg
    end
    if  ~isempty(rips) & ...
            isfield(SleepState.ints,'NREMstate') & ...
            ~isempty(ls)
        nCells_ls = length(ls.times);
        [sleep_ls] = double(InIntervals(rips.peaks,SleepState.ints.NREMstate));
        state_ls = [state_ls;sleep_ls];
        for r = 1:length(rips.timestamps)
            spks = (InIntervals(ls.spindices(:,1),[rips.timestamps(r,1)-.025 rips.timestamps(r,2)+.025]));
            ls_nSpks(r) = length(unique(ls.spindices(spks,2)));
        end
        ls_prop = [ls_prop,ls_nSpks./nCells_ls]; clear ls_nSpks
    end
%     cd /home/david/datasets/ripples_LS 
    cd E:\datasets\ripples_LS



r = unique(region);

for rr = 1:length(r)
   subplot(2,2,rr)
   
   s = find(state==1);
   histogram(log(proportion(intersect(s,find(region==r(rr))))),[-5:.25:1],'Normalization','pdf','FaceColor','k')
   s = find(state==0);
   hold on
   histogram(log(proportion(intersect(s,find(region==r(rr))))),[-5:.25:1],'Normalization','pdf','FaceColor','b')
   title(r(rr))
   
end

   subplot(2,2,4)
   s = find(state_ls==1);
   histogram(log(ls_prop(s)),[-5:.25:1],'Normalization','pdf','FaceColor','k')
   hold on
   s = find(state_ls==0);
   histogram(log(ls_prop(s)),[-5:.25:1],'Normalization','pdf','FaceColor','b')
   hold off
   pause(.1)
clf
end