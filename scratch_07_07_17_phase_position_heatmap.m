
% [firingMaps] = bz_firingMap1D(spikes.times,behavior,lfp,4);
d = dir('*firingMaps*');
spikes = bz_GetSpikes;
load(d.name);
d = dir('*.behavior.*');
load(d.name);
for cell =1:length(spikes.times)
    if strcmp(spikes.region{cell},'ls')
    for condition = 1:length(firingMaps.rateMaps)
        times = [];
        f = find(behavior.events.trialConditions==condition);
        for t=1:length(f)
        times = [times;Restrict(spikes.times{cell},behavior.events.trialIntervals(f(t),:))];
        end
        p = firingMaps.phaseMaps{condition}{cell};
       
    if ~isempty(p) & length(p)>30
%     col = (cos([p(:,end)']'));
    col = p(:,end-1);
    col_orig = col;
    [a ind] = sort(p(:,1));
    col = Smooth(col(ind),length(p)./10);
    %         col = wrap(circ_smoothTS(p(ind,end),round(length(p)./10),'method','median'));
    %         col_sorted = col(ind);
    %         col = minmax_norm([-1 1 col']);
    %         col = col(3:end);
    col_rate = sort_back(col,ind);
    
%     col = (sin([p(:,end)']'));
    [a ind] = sort(p(:,1));
    col = wrap(circ_smoothTS(p(ind,end),round(length(p)./10),'method','median'));
% % col = wrap(p(ind,end));
    col_phase = sort_back(col,ind);
%     t=trials{condition};
%     tt=[];
%     for i=1:length(trials)
%     tt=[tt;t{i}];
%     end
%     for i = length(keep)
%     f = find(p(:,5)<times(i));
    f = 1:length(times);
    %             if indices(i) < 12 || i == length(keep)
%     figure
    %             jetwrap = vertcat(jet,flipud(jet));
    %             colormap(jetwrap)
    %         set(gcf,'Visible', 'off'); 
%     colormap(parula)
%     imshow(keep{i})
    %     scatter(tt(1:i,2),tt(1:i,1),'.k')
%     axis([120 510 50 450 ])
%     hold on
% 
%     for jj = 1:length(f)
%     %        [a b] = min(abs(pos(:,1)-p(f(j),5)));
%     %        scatter(mean(pos(b,[2 4])'),mean(pos(b,[1 3])'),'.r')
%     scatter(p(f(jj),3),p(f(jj),4),[],col_rate(f(jj)),'.')
%     end
%     caxis([0 20])
    clf
    jetwrap = vertcat(jet,flipud(jet));
    colormap(jetwrap)
    axis square
    %         set(gcf,'Visible', 'off'); 
    %             imshow(keep{i})
    %     scatter(tt(1:i,2),tt(1:i,1),'.k')
%     axis([-1000 1000 -1000 1000])
    for i =1:length(behavior.events.trialConditions)
       behavior.events.trials{i}.x = makeLength(behavior.events.trials{i}.x,250);
       behavior.events.trials{i}.y = makeLength(behavior.events.trials{i}.y,250);
       behavior.events.trials{i}.z = makeLength(behavior.events.trials{i}.z,250);
    end
%     bz_plotTrials(behavior,'condition',condition,'color',[.5 .5 .5],'endPoints',false)
    hold on
%     for jj = 1:length(f)
%     %        [a b] = min(abs(pos(:,1)-p(f(j),5)));
%     %        scatter(mean(pos(b,[2 4])'),mean(pos(b,[1 3])'),'.r')
%     scatter(p(f(jj),3),p(f(jj),4),[],col_phase(f(jj)),'.')
%     end
    scatter(p(:,3),p(:,4),[],col_phase,'.')
    caxis([-pi pi])
    if isempty(f)
    text(130,60,['trials: 1'],'color','w')
    else
    text(130,60,['trials: ' num2str(sum(behavior.events.trialConditions==condition))],'color','w')
    end
    set(gca,'color','black')
    colorbar
    xlabel([ spikes.sessionName ' cell #: ' num2str(cell)])
    
% pause
    hold off
%     end
pause
try
%      saveFigure(['/home/david/Dropbox/tmp/' ...
%             spikes.sessionName '_' num2str(cell) '_' num2str(condition) '.fig'])
close
catch
    close all
end
pause(.1)
    clf
    end
    end
    end
end