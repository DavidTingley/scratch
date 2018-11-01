d = dir('*201*');
c=1;
region = [];
proportion =[];
state = [];
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
    
    if ~isempty(hpc) & ~isempty(rips) & isfield(SleepState.ints,'NREMstate') & ~isempty(ls)
        
        nCells = length(hpc.times);
        nCells_ls = length(ls.times);
        
        [sleep] = double(InIntervals(rips.peaks,SleepState.ints.NREMstate));

        for r = 1:length(rips.timestamps)
            spks = (InIntervals(hpc.spindices(:,1),[rips.timestamps(r,1)-.025 rips.timestamps(r,2)+.025]));
            nSpks(r) = length(unique(hpc.spindices(spks,2)));
            reg(r) = sum(double(hpc.region{1}));
            ls_nSpks(r) = length(Restrict(ls.spindices(:,1),rips.timestamps(r,:)));
        end


        state = [state;sleep];
        region = [region, reg];
        proportion = [proportion, nSpks./nCells]; clear nSpks reg
        ls_prop = [ls_prop,ls_nSpks./nCells_ls]; clear ls_nSpkes
        
    end
    cd /home/david/datasets/ripples_LS 
end



r = unique(region);

for rr = 1:length(r)
   subplot(2,2,rr)
   
   s = find(state==1);
   histogram(log(proportion(intersect(s,find(region==r(rr))))),[-5:.5:1],'FaceColor','k')
   s = find(state==0);
   hold on
   histogram(log(proportion(intersect(s,find(region==r(rr))))),[-5:.5:1],'FaceColor','b')
   title(r(rr))
   
end
%    