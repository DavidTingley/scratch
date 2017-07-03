% this script tries to examine the information content of phase and rate
% codes in HPC at different temporal/spatial integration windows.  by smoothing
% with different spatial scales, we can examine how information is
% gained/lost in these two codes. 
for cond = 1:8
    c=countMap{cond};
    phases_trial = zeros(length(phaseMap{cond}),size(rateMap{cond},2),(size(rateMap{1},3)));

    for j=1:size(rateMap{cond},2)
            for cell = 1:length(phaseMap{cond})
            for k = 1: size(rateMap{1},3)
                if ~isempty(phaseMap{cond}{cell})
                f = find(phaseMap{cond}{cell}(:,1)==k); 
                ff = find(phaseMap{cond}{cell}(:,2)==j);
                if ~isempty(intersect(f,ff))
                    phases_trial(cell,j,k) = circ_mean(phaseMap{cond}{cell}(intersect(f,ff),end));
                end
                end
            end
            end
    end
    
for i=1:5
    for j=1:length(spikes.times)
        for k=1:size(rateMap{cond},2)
            phases_smooth(j,k,:)=circ_smoothTS(squeeze(phases_trial(j,k,:)),i,'method','mean','exclude',0);
        end
    end
    phases_smooth(phases_smooth==0)=nan;
    [phases_smooth, EDGES] = discretize(phases_smooth,-pi:.1:pi);
    phases_smooth(isnan(phases_smooth))=0;
    
    [track_info_phase(i,:,:) pos_info] = Info_Analysis(phases_smooth,5,0);

    [track_info_rate(i,:,:) pos_info] = Info_Analysis(c,2,i);





    subplot(2,2,1);
    imagesc(squeeze(mean(track_info_rate,2)));
    subplot(2,2,2);
    imagesc(squeeze(mean(track_info_phase,2)));

    subplot(2,2,3)
    plot(mean(squeeze(track_info_phase(:,80,:)),2),'g');
    hold on
    plot(mean(squeeze(track_info_rate(:,80,:)),2),'r');
    hold off

    subplot(2,2,4)
    f = find(spikes.shankID<5);
    ff = find(spikes.shankID>4);
    
    plot(squeeze(nanmean(nanmean(track_info_rate(:,f,:),3),2)),'.r')
    hold on
    plot(squeeze(nanmean(nanmean(track_info_phase(:,f,:),3),2)),'.g')
    plot(squeeze(nanmean(nanmean(track_info_rate(:,ff,:),3),2)),'r')
    plot(squeeze(nanmean(nanmean(track_info_phase(:,ff,:),3),2)),'g')
    hold off

%         subplot(2,2,2)
%         imagesc(squeeze(phases_smooth(80,:,:)))
    i
    pause(.1)
end
    info_rate{cond} = track_info_rate;
    info_phase{cond} = track_info_phase;
end

