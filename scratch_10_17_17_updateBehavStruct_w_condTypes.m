d  = dir('*201*'); 

for ii=1:length(d)
    cd(d(ii).name)
    load([d(ii).name '.behavior.mat'],'behavior')
    if ~isfield(behavior.events,'conditionType')
    if strcmp(behavior.description ,'central alternation')
        for cond = 1:length(behavior.events.map)
        behavior.events.conditionType{cond} = 'central';
        end
    elseif strcmp(behavior.description ,'wheel alternation')
        for cond = 1:length(behavior.events.map)
        behavior.events.conditionType{cond} = 'wheel';
        end
    elseif strcmp(behavior.description ,'linear')
        for cond = 1:length(behavior.events.map)
        behavior.events.conditionType{cond} = 'linear';
        end
    else
        bz_plotTrials(behavior)
        for cond = 1:length(behavior.events.map)
           s(cond) = input(['trial type for ' num2str(cond) ': '],'s');
           if strcmp(s(cond),'w')
               behavior.events.conditionType{cond} = 'wheel';
           end
           if strcmp(s(cond),'c')
               behavior.events.conditionType{cond} = 'central';
           end
           if strcmp(s(cond),'l')
               behavior.events.conditionType{cond} = 'linear';
           end
        end 
    end

    save([d(ii).name '.behavior.mat'],'behavior')
    end
    clear s behavior
    clf
    
    cd /home/david/datasets/lsDataset
end
