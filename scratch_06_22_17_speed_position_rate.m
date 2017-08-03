% function [] = scratch6_22_17_speed_position_rate()
% here we are examing the relationship between running speed, place fields
% and their instantaneous rates.  preliminary analysis suggests that the
% correlation between firing rate and running speed VARIES as a function of
% distance the animals goal location...



    
% % to be run from within a recording folder
d= dir('*');
for rec = 3:length(d)
    cd(d(rec).name)
% try
    clearvars -except rec d rate_vel rate_acc center_of_mass  rate_vel_pval beh UID spk ratemaps veloc ratemaps fast_trials slow_trials
    c=1;
    % get the data for a recording
    xml = LoadParameters;cd 
    if exist([xml.FileName '.spikes.cellinfo.mat'])
        load([xml.FileName '.spikes.cellinfo.mat'])
    else
        spikes = bz_GetSpikes('region','hpc');    
    end
        
    load([xml.FileName '.behavior.mat'])
    load([xml.FileName '.sessionInfo.mat'])
    
    lfp = bz_GetLFP(sessionInfo.thetaChans(end));%,'intervals',[behavior.timestamps(1) behavior.timestamps(end)]);
                                                 % I need to fix this and
                                                 % other functions *1250
                                                 % breaks when only an
                                                 % interval is loaded in
                                                 % for lfp struct
    if strcmp(behavior.trackingType,'optitrack')
        tau = 5;
        maxFieldWidth = 100;
        minFieldWidth = 8;
    else
       tau = 2; 
       maxFieldWidth = 40;
       minFieldWidth = 3;
    end
    [rateMap countMap opk_rate_vel_corruMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,tau);
    
    figure(1)
    hpc = [];
    for i=1:length(spikes.times)
       if strcmp(spikes.region{i},'hpc')
           hpc(i) = 1;
       else
           hpc(i) = 0;
       end
    end
    % find place fields
    for i =1:length(unique(behavior.events.trialConditions))
        fields{i} = bz_getPlaceFields1D(rateMap{i},'minPeakRate',5,'maxFieldWidth',maxFieldWidth,'minFieldWidth',minFieldWidth);
    end

    %% END LOADING DATA
    % pk_rate_vel_corrumulate stats on fields and behavior
    for i=1:length(unique(behavior.events.trialConditions))
        f = find(behavior.events.trialConditions==i);
    for j=find(hpc)
    if ~isempty(fields{i}{j}) & length(f) > 20
        
        for ff = 1:length(fields{i}{j})
        for ii = 1:length(find(behavior.events.trialConditions==i))
            for v = 1:length(behavior.events.trials{f(ii)}.x)-1
                vel(v) = sum(sqrt(([behavior.events.trials{f(ii)}.x(v) behavior.events.trials{f(ii)}.y(v)] ...
                         - [behavior.events.trials{f(ii)}.x(v+1) behavior.events.trials{f(ii)}.y(v+1)]).^2));
            end
            vel = smoothts(squeeze(vel),'e',maxFieldWidth); vel(end+1) = vel(end);
            acc = (squeeze(diff(vel))); acc(end+1) = acc(end);
            for v=1:length(behavior.events.map{1}.x)
                vv(ii,v) = mean(vel(find(behavior.events.trials{f(ii)}.mapping==v)));
                aa(ii,v) = mean(acc(find(behavior.events.trials{f(ii)}.mapping==v)));
            end;
            [peakR(ii) loc] = max(rateMap{i}(j,ii,fields{i}{j}{ff}.start:fields{i}{j}{ff}.stop));
            loc = loc + fields{i}{j}{ff}.start;
%             [peakR(ii)] = max(rateMap{i}(j,ii,fields{i}{j}{ff}.COM));
%             if loc > size(rateMap{1},3)-5
%                 loc = size(rateMap{1},3)-5;
%             elseif loc < 5
%                 loc = 5;
%             end
%             peakR(ii) = mean(rateMap{i}(j,ii,loc-4:loc+4)); 
            num_spks{i}{j}(ii) = sum(countMap{i}(j,ii,:));
            avg_vel{i}{j}(ii) = mean(vel);
            velat(ii) =vv(ii,fields{i}{j}{ff}.COM);
            accat(ii) = aa(ii,fields{i}{j}{ff}.COM);
            clear vel acc aa vv;
        end
        peakR(peakR==0)=nan;  % for now we remove 0 rate trials
        peak_rates{i}{j}(ff,:) = peakR;
        velocities{i}{j}(ff,:) = velat; clear peakR velat;
        accelerations{i}{j}(ff,:) = accat; clear accat;
        peak_loc{i}{j}(ff) = fields{i}{j}{ff}.peakLoc; 
        chan{i}{j}(ff) = spikes.maxWaveformCh(j);
        peak_com{i}{j}(ff) = fields{i}{j}{ff}.COM; 
        shank{i}{j}(ff) = spikes.shankID(j);
        end
    else
        peak_rates{i}{j} = nan;
        velocities{i}{j} = nan; clear peakR velat;
        peak_loc{i}{j} = nan; 
        chan{i}{j} = nan;
        peak_com{i}{j} = nan; 
        shank{i}{j} = nan;
        num_spks{i}{j}=nan;
        avg_vel{i}{j}=nan;
    end
    end
    end
    
    % below can be used to plot single cells or accumulate across
    % behavioral conditions and place fields

    for i=1:length(unique(behavior.events.trialConditions))
    for j=1:length(peak_rates{i})
    for ff = 1:length(fields{i}{j})
        if ~isempty(peak_rates{i}{j}) & ~isempty(fields{i}{j}) & length(num_spks{i}{j})>1
        subplot(2,2,1)
        scatter(peak_rates{i}{j}(ff,:),velocities{i}{j}(ff,:),'.k')
        subplot(2,2,2)
        histogram(peak_rates{i}{j}(ff,:)./velocities{i}{j}(ff,:),20)
        subplot(2,2,3)

        [pk_rate_vel_corr(c) pk_rate_vel_pval(c)] = corr(peak_rates{i}{j}(ff,:)',velocities{i}{j}(ff,:)','rows','complete');
        [num_spk_vel_corr(c) num_spk_vel_pval(c)] = corr(num_spks{i}{j}',velocities{i}{j}(ff,:)','rows','complete');
        [pk_rate_acc_corr(c) pk_rate_acc_pval(c)] = corr(peak_rates{i}{j}(ff,:)',accelerations{i}{j}(ff,:)','rows','complete');
        [num_spk_acc_corr(c) num_spk_acc_pval(c)] = corr(num_spks{i}{j}',accelerations{i}{j}(ff,:)','rows','complete');
        
%         for iter = 1:100
%             [pk_rate_vel_corr(c) pk_rate_vel_pval(c)] = corr(peak_rates{i}{j}(ff,:)',velocities{i}{j}(ff,:)','rows','complete');
%         end

        shankID(c) = shank{i}{j}(ff);
        chanID(c) = chan{i}{j}(ff);
        pk_loc(c) = peak_loc{i}{j}(ff);
        com(c) = peak_com{i}{j}(ff);  % why isn't this being returned? later... 
        condition(c) = i;
        
        cellUID(c) = j;
%       title(pk_rate_vel_corr(c)); 
%       [a b] = sort(velocities{i}{j}(ff,:));
%       imagesc(squeeze(rateMap{i}(j,b,:)))
%       subplot(2,2,4)
%       plot(spikes.rawWaveform{j})
        [a b] = sort(velocities{i}{j}(ff,:));
        slow = 1:ceil(prctile(b,50));
        %             med = slow(end)+1:ceil(prctile(b,66));
        fast = slow(end)+1:length(b);
        plot(squeeze(mean(rateMap{i}(j,slow,:))),'g')
        hold on
        %             plot(squeeze(mean(rateMap{i}(j,med,:))),'k')
        plot(squeeze(mean(rateMap{i}(j,fast,:))),'r'); title(j);

        fas_trials(c,:)= squeeze(mean(rateMap{i}(j,fast,:))); 
        slo_trials(c,:)= squeeze(mean(rateMap{i}(j,slow,:)));
            pause
        else
        pk_rate_vel_corr(c) = nan;
        pk_rate_vel_pval(c) = nan;
        num_spk_vel_corr(c) = nan;
        num_spk_vel_pval(c) = nan;
        pk_rate_acc_corr(c) = nan;
        pk_rate_acc_pval(c) = nan;
        num_spk_acc_corr(c) = nan;
        num_spk_acc_pval(c) = nan;
        shankID(c) = nan;
        chanID(c) = nan;
        cellUID(c) = nan;
        pk_loc(c) = nan;
        com(c) = nan;  % why isn't this being returned? later... 
        condition(c) = i;
        end
    c=c+1;
    end
    end
    end
    
    % for visualization
%     for j=find(hpc)
%         for i=1:length(unique(behavior.events.trialConditions))
%             for ff = 1:length(fields{i}{j})
%             if ~isempty(velocities{i}{j}) & ~isnan(velocities{i}{j})
%                 if ~isempty(fields{i}{j})
%             subplot(5,2,10)
%             [a b] = sort(velocities{i}{j}(ff,:));
%             slow = 1:ceil(prctile(b,50));
% %             med = slow(end)+1:ceil(prctile(b,66));
%             fast = slow(end)+1:length(b);
%             plot(squeeze(mean(rateMap{i}(j,slow,:))),'g')
%             hold on
% %             plot(squeeze(mean(rateMap{i}(j,med,:))),'k')
%             plot(squeeze(mean(rateMap{i}(j,fast,:))),'r'); title(j);
%             
%             f(c,:)= squeeze(mean(rateMap{i}(j,fast,:))); 
%             s(c,:)= squeeze(mean(rateMap{i}(j,slow,:)));
%             
%             subplot(5,2,i)
%             imagesc(squeeze(rateMap{i}(j,b,:)))
%             if ~isempty(peak_rates{i}{j})
%                [a p]=  corr(peak_rates{i}{j}(ff,:)',velocities{i}{j}(ff,:)','rows','complete');
%             title([a p])
%             end
%                 end
%             end
%             
%             end
%         end
% %         pause; clf
%     end
    
%     clear peak_rates velocities
%     title(sessionInfo.FileName)
    figure(101)
    subplot(2,2,1)
    scatter(Bin(com,[0 size(rateMap{1},3)],80),pk_rate_vel_corr,'.')
    ylabel('pk rate / veloc correlation')
    xlabel('relative position of field (0-start, 80-end)')
    hold on
    title('COM vs corr, VELOCITY')
    axis([0 80 -1 1])
    
    subplot(2,2,2)
    scatter(Bin(com,[0 size(rateMap{1},3)],80),num_spk_vel_corr,'.')
    ylabel('num spikes / veloc correlation')
    xlabel('relative position of field (0-start, 80-end)')
    hold on
    title('com vs corr, VELOCITY')
    axis([0 80 -1 1])
    
    subplot(2,2,3)
    scatter(Bin(com,[0 size(rateMap{1},3)],80),pk_rate_acc_corr,'.')
    ylabel('peak rate / veloc correlation')
    xlabel('relative position of field (0-start, 80-end)')
    hold on
    title('COM vs peak rate corr, ACCELERATION')
    axis([0 80 -1 1])
    
    subplot(2,2,4)
    scatter(Bin(com,[0 size(rateMap{1},3)],80),num_spk_acc_corr,'.')
    ylabel('num spks / veloc correlation')
    xlabel('relative position of field (0-start, 80-end)')
    hold on
    title('COM vs num spks corr, ACCELERATION')
    axis([0 80 -1 1])
    
%     rate_vel{rec} = pk_rate_vel_corr;
%     rate_acc{rec} = pk_rate_acc_corr;
%     center_of_mass{rec} = Bin(com,[0 size(rateMap{1},3)],80);
%     rate_vel_pval{rec} = pk_rate_vel_pval;
%     beh{rec} = behavior;
%     UID{rec}= cellUID;
%     spk{rec} = spikes;
%     veloc{rec} = velocities;
%     ratemaps{rec} = rateMap;
%     fast_trials{rec}=fas_trials;
%     slow_trials{rec}=slo_trials;
    
%     figure(102)
%     subplot(2,2,1)
    f = (fast_trials{i}');
    s = (slow_trials{i}');
    for i =1:size(f,1)
       ff(i,:)=minmax_norm(f(i,:)); 
       ss(i,:)=minmax_norm(s(i,:)); 
       di(i,:) = zscore(f(i,:)-s(i,:));
    end
    [a b o]= sort_cells(ff,ss,1);
    imagesc(ff(o,:));
    
    subplot(2,2,4)
    imagesc(di(o,:))
    
    pause(.2)
    cd /home/david/datasets/speedRate
end
save('speed_acc_metadata.mat')



% catch
% end
% end