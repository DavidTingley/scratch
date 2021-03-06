[b a] = butter(4,[120/625 200/625],'bandpass');
d = dir('*');
for r = 3:length(d)
    if d(r).isdir
        cd(d(r).name)
        if exist([d(r).name '.ripplesALL.event.mat']) & ~isempty(dir('*pos')) & ~isempty(dir('*clu*'))
            figure
        load([d(r).name '.ripplesALL.event.mat'])
        lfp = bz_GetLFP([27 ripples.rippleChan]);
        spikes = bz_GetSpikes('savemat',true,'forcereload',true,'getwaveforms',true);
        beh = dir('*pos');
        xml = LoadParameters;
        
        pos = importdata(beh.name);
%       if behav do this
        if exist('maps.mat') | ~isempty(dir('behav*.mat'))
            try
                load('maps.mat')
                for i = 1:length(rateMap)
                [track_info{i} pos_info{i}] = Info_Analysis(rateMap{i}*1000,2);
                end
                info = mean(cell2mat(track_info)');
            catch
                t = dir('behav*mat');
                load(t(1).name)
                if length(trials) < 12
                    if size(pos,2) == 5
                        [rateMap occuMap phaseMap] = spaceRateMap_old(spikes.times,pos,map,mapping,trials,[lfp.timestamps ,double(lfp.data(:,1))]);
                    else
                        [rateMap occuMap phaseMap] = spaceRateMap(spikes.times,pos,map,mapping,trials,lfp);
                    end
                for i = 1:length(rateMap)
                [track_info{i} pos_info{i}] = Info_Analysis(rateMap{i}*1000,2);
                end
                info = mean(cell2mat(track_info)');
                else
                        info = nan(length(spikes.times),1);
                end
            end
        


%         else
    %       else do this
%             for i = 1:length(spikes.times)
%                [fr_map{i} stats{i}] = FiringMap(pos(:,[1 8 10]),spikes.times{i},'nBins',200); 
%                info(i) = stats{i}.specificity;
%             end 
%         end
        pow(:,1) = FiltFiltM(b,a,double(lfp.data(:,1)));
        pow(:,2) = FiltFiltM(b,a,double(lfp.data(:,2)));
        
       
        
        
        
        
        for i=1:length(ripples.times)
        p(i) = max(pow(round(ripples.peaks(i)*1250)-40:round(ripples.peaks(i)*1250)+40,1));
        pp(i) = max(pow(round(ripples.peaks(i)*1250)-40:round(ripples.peaks(i)*1250)+40,2));
        end
        subplot(2,2,1)
        scatter(p,pp,'.')
        hold on
        f = find(p>mean(p)+std(p)*2);
        ff = find(p<=mean(p)+std(p)*2);
% for i=1:length(ripples.times)
% p(i) = max(pow(round(ripples.peaks(i)*1250)-40:round(ripples.peaks(i)*1250)+40,1));
% pp(i) = max(pow(round(ripples.peaks(i)*1250)-40:round(ripples.peaks(i)*1250)+40,2));
% end
        
        
        for i=1:length(ripples.times)
        for j=1:length(spikes.times)
        sp = find(spikes.times{j}>ripples.peaks(i)-.04);
        sp2 = find(spikes.times{j}<ripples.peaks(i)+.04);
        rate(i,j,ceil((spikes.times{j}(intersect(sp,sp2))-ripples.peaks(i))*1250)+50)=1;
        end
        end
        delta_rate = (squeeze(mean(mean(rate(f,:,:),3))-mean(mean(rate(ff,:,:),3)))) ./ squeeze(mean(mean(rate,3)));

        ls_ind=[]; hpc = [];
        for i=1:length(spikes.times)
            if strcmp(spikes.region{i},'ls')
                ls_ind = [ls_ind i];
            elseif strcmp(spikes.region{i},'hpc')
                hpc = [hpc i];
            end
        end
        
        subplot(2,2,2)
        if length(hpc) > 0
        plot(info(hpc),delta_rate(hpc),'.k');
        ylabel('coupled - decoupled')
        xlabel('spatial info')
        hold on
        subplot(2,2,3)
        negsemilogx(delta_rate(hpc),info(hpc))
        hold on
        end
        subplot(2,2,4)
        if length(ls_ind) > 0
        plot(info(ls_ind),delta_rate(ls_ind),'.r');
        hold on
        end
        
        end
        clear rate p pp pow info delta* ls_ind hpc track_info pos_info *Map
%         savefig([xml.FileName(1:4) '.fig'])
        end
        cd ..
        pause(.5)
        
    end
end

