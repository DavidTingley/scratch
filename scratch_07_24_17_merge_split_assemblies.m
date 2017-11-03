d  = dir('*201*');
% % ii=38;
% 
for ii=1:length(d)
    for cond = 1:20
        if exist(['/mnt/packrat/userdirs/david/zpool1/results/' d(ii).name '_condition_' num2str(cond) '.mat'])
            dat = load(['/mnt/packrat/userdirs/david/zpool1/results/' d(ii).name '_condition_' num2str(cond) '.mat']);
            dev{cond} = dat.dev{cond};
            devControl{cond} = dat.devControl{cond};
            pairs = dat.pairs;
            velocities = dat.velocities;
            phasetrains = dat.phasetrains;
            spktrains = dat.spktrains;
            coords = dat.coords;
            if cond ~= 1
            if ~isempty(dev{cond}) & isempty(dev{cond-1})
               warning(['missing ' d(ii).name ' ' num2str(cond)]) 
            end
            end
        end   
    end
    if exist('dat')
        save(['/home/david/datasets/lsDataset/' d(ii).name '/assembliesCrossRegion_split_w_theta_' date '.mat']);
        clear dev* dat* *trains velocities coords
    end
    ii
end

% d  = dir('*201*');
% % % ii=38;
% % 
% for ii=1:length(d)
%     for cond = 1:20
%         if exist(['/mnt/packrat/userdirs/david/zpool1/results/' d(ii).name '_w_positionFits_condition_' num2str(cond) '.mat'])
%             load(['/mnt/packrat/userdirs/david/zpool1/results/' d(ii).name '_w_positionFits_condition_' num2str(cond) '.mat']);
%             for res = 1:size(peerPredictionFits.results,1)
%                 
%                 peerPredictionFits.results.fits.yfit              
% 
%             end    
%         
%         end   
%     end
%     if exist('dat')
%         save(['/home/david/datasets/lsDataset/' d(ii).name '/assembliesCrossRegion_split_w_theta_w_positionFits_.mat']);
%         clear dev* dat* *trains velocities coords
%     end
% end