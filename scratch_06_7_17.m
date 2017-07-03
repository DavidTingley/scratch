d = dir('*');
for r = 3:length(d)
    if d(r).isdir
        cd(d(r).name)
        if ~exist([d(r).name '.pos'])  & ~isempty(dir('Session*'))
        Process_ConvertOptitrack2Pos(d(r).name)
        b = dir('*pos');
        pos = importdata(b.name);
        subplot(6,6,r)
        scatter(pos(:,8),pos(:,10));
        end
        cd ..
    end
end
