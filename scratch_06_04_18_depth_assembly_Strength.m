d = dir('*201*');
% figure
count = 1;
clf
for i=1:length(d)
cd(d(i).name)
if exist('assembliesCrossRegion_split_w_theta.mat')
load('assembliesCrossRegion_split_w_theta.mat','pairs','dev*')
ripples = bz_LoadEvents(pwd,'CA1Ripples');
spikes = bz_GetSpikes('noprompt',true);
if ~isempty(ripples)
    for spk = 1:length(spikes.times)
       spikes.times{spk} = Restrict(spikes.times{spk},[ripples.peaks-.25 ripples.peaks+.25]); 
    end
    [times groups] = spikes2sorted(spikes.times);
    [ccg t] = CCG(times,groups,'binSize',.001);
    ls_idx = find(strcmp(spikes.region,'ls'));

    for cell=1:length(spikes.times)
       if strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca1') 
%         f = find(pairs(:,2)==cell);
%         imp = nan(length(dev),length(pairs));
%         for cond = 1:length(dev)
%             for pair = 1:size(dev{cond},2)
%                  [a b] =  min(dev{cond}(:,pair));
%                  [aa bb] = min(mean(devControl{cond}(:,pair,:),3));
%                  imp(cond,pair) = (a-mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(:,pair,:),3)));
%             end
%         end
        for l = 1:length(ls_idx)
            cellCCG(count,:) = zscore(fastrms(squeeze(mean(ccg(:,cell,ls_idx(l)),3)),7));
%         if ~isempty(f)
%             cellContrib(count) = nanmean(nanmean(imp(:,f)));
            cellDepth(count) = spikes.chanDepthRelative_CA1PYR_wav(cell);
            count = 1+count;
%         else
%             cellContrib(count) = nan;
%             cellDepth(count) = nan;
%             count = 1+count;
        end
%         clear imp
%        end
%        else
%            cellContrib(count) = nan;
%            cellDepth(count) = nan;
%            count = 1+count;
       end
    end
    if isfield(spikes,'chanDepthRelative_CA1PYR_wav')
%         subplot(2,2,1)
%         plot(spikes.chanDepthRelative_CA1PYR_wav,cellContrib,'.k')
%         hold on
        subplot(2,2,2)
        histogram(cellDepth,-100:5:60)
        subplot(2,2,3)
        for ii = -100:60
           idx = find(cellDepth>=ii-5 & cellDepth<=ii+5); 
           if ~isempty(idx)
%                plot(ii,nanmean(cellContrib(idx)),'.r')
%                hold on
               cc(ii+101,:) = (nanmean(cellCCG(idx,:)));
           else
               cc(ii+101,1:2001) = nan;
           end
        end
        for ii=1:161
        plot(cc(ii,:)*5+ii/2,'k')
        hold on
        end
        axis([750 1250 0 81])
        hold off
        subplot(2,2,4)
        imagesc(-100:60,-1000:1000,cc')
        axis([-100 60 -250 250])
        subplot(2,2,1)
        plot(-100:60,nanmean(cc(:,1002:1015),2)-nanmean(cc(:,987:1000),2))
        pause(.1)
    end
%     clear cellContrib
end
end

cd /home/david/datasets/ripples_LS/
end