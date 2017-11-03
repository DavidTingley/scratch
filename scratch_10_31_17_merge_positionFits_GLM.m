d  = dir('*201*');
% % ii=38;
% 
count = 1;
for ii=1:length(d)
    cd(d(ii).name)
    for cond = 1:20
        if exist(['/mnt/packrat/userdirs/david/zpool1/results/' d(ii).name '_w_positionFits_condition_' num2str(cond) '.mat'])
            load(['/mnt/packrat/userdirs/david/zpool1/results/' d(ii).name '_w_positionFits_condition_' num2str(cond) '.mat']);
            for res = 1:size(peerPredictionFits.results,1)
                for pos = 1:201
                    f = find(abs(peerPredictionFits.results.fits(res).position-pos)<2);
                    se(count,pos) = single(mean((peerPredictionFits.results.fits(res).yfit(f)-peerPredictionFits.results.fits(res).rate(f)).^2));
                end
                if isempty(sessionInfo.ca3)
                    region(count) = 1;
                else
                    region(count) = 3;
                end
                
                count=1+count;
            end    
        end   
        plot(mean(se))
        pause(.1)
    end
    cd /home/david/datasets/lsDataset
%     if exist('dat')
%         save(['/home/david/datasets/lsDataset/' d(ii).name '/assembliesCrossRegion_split_w_theta_w_positionFits_.mat']);
%         clear dev* dat* *trains velocities coords
%     end
end