d  = dir('*201*');
% % ii=38;
% 
for ii=1:length(d)
    for cond = 1:20
        if exist(['/home/david/datasets/results/' d(ii).name '_condition_' num2str(cond) '.mat'])
            dat = load(['/home/david/datasets/results/' d(ii).name '_condition_' num2str(cond) '.mat']);
            dev{cond} = dat.dev;
            devControl{cond} = dat.devControl;
            pairs = dat.pairs;
            velocities = dat.velocities;
            phasetrains = dat.phasetrains;
            spktrains = dat.spktrains;
            coords = dat.coords;
        end   
    end
    if exist('dat')
        save(['/home/david/datasets/lsDataset/' d(ii).name '/assembliesCrossRegion_split_w_theta.mat']);
        clear dev* dat* *trains velocities coords
    end
end