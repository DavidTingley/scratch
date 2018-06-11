d = dir('*201*');
% figure
clf
for i=1:95
cd(d(i).name)
if exist('assembliesCrossRegion_split_w_theta.mat')
load('assembliesCrossRegion_split_w_theta.mat','pairs','dev*')

spikes = bz_GetSpikes('noprompt',true);

    for cell=1:length(spikes.times)
       if strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca1') 
        f = find(pairs(:,2)==cell);
        imp = nan(length(dev),length(pairs));
        for cond = 1:length(dev)
            for pair = 1:size(dev{cond},2)
                 [a b] =  min(dev{cond}(:,pair));
                 [aa bb] = min(mean(devControl{cond}(:,pair,:),3));
                 imp(cond,pair) = (a-mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(:,pair,:),3)));
            end
        end
        if ~isempty(f)
            cellContrib(cell) = nanmean(nanmean(imp(:,f)));
        else
            cellContrib(cell) = nan;
        end
        clear imp
       else
           cellContrib(cell) = nan;
       end
    end
    if isfield(spikes,'chanDepthRelative_CA1PYR_wav')
        plot(spikes.chanDepthRelative_CA1PYR_wav,cellContrib,'.k')
        hold on
        pause(.1)
    end
    clear cellContrib
end

cd /home/david/datasets/ripples_LS/
end