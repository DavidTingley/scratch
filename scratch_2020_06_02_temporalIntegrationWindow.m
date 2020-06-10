files = dir('_*mat');
iter = 1;
for f =1:length(files)
    load(files(f).name,'count','at','rt','ripples','isig_levels','absTime','idx')
    
    idx = makelength(idx,length(idx)*5*2);
    idx = round(idx .* 5*2);
    
    isa = [0,diff(isig_levels)];
    for i=1:length(isa)
        isa_1min(i*5*2:i*5*2+5*2)=isa(i);
    end
    abso = (makeLength(absTime,length(isa_1min)));
    
    c=1;
    for ii=1:length(ripples)
        if ~isempty(ripples{ii})
        for i=1:length(ripples{ii}.peaks)
            [blah b]= min(abs(rt{ii}-ripples{ii}.peaks(i)));
            a(c) = at{ii}(b);
            c=1+c;
        end
        end
    end
    absolute = a; clear a
    
    parfor i=1:length(abso)
        ind = find(InIntervals(absolute,[abso(i)-1.15741277113557e-05*30 abso(i)]));
        count(i) = length(ind);
    end
    counts{f} = count;
    idxs{f} = idx;
    isa_1mins{f} = isa_1min;
    
    for i=1:50
        cc(f,i,:) = (ccgBinned(movsum(count(idx),i),isa_1min(idx),200));
    end
    clear count isa_1min
end
