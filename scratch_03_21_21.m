

d = dir('*201*'); % get a list of recordings to run through (N = 98 in LS dataset, ripples_LS 2020 paper dataset)

for i=1:length(d)
   cd(d(i).name)
   sessionInfo = bz_getSessionInfo;
   
    if exist([d(i).name '.CA1Ripples.events.mat']) & exist([d(i).name '.spikes.cellinfo.mat']) % if we have spikes and ripples for this rec
       
       load([d(i).name '.CA1Ripples.events.mat'])
       spikes = bz_GetSpikes('region','hpc');
       if isempty(spikes)
           spikes = bz_GetSpikes('region','ca1'); % I was lazy and labelled diff recordings CA1 or HPC.. and still haven't changed things
       end
       if ~isempty(spikes) % @Dan, we should have fixed this in bz_GetSpikes ages ago to return an empty struct, rather than just empty..
       if length(spikes.times) > 10 % threshold, need at least 10 cells to use recording
           
           rate = zeros(ceil(spikes.spindices(end,1)*1000),1); % 1ms bins, for the whole recording
           rateTS = 0:1/1000:spikes.spindices(end,1);

           for spk = 1:length(spikes.spindices)
                rate(ceil(spikes.spindices(spk,1)*1000)) = rate(ceil(spikes.spindices(spk,1)*1000)) + 1; % really a count, but it doesnt matter
           end
           populationZscoredRate = zscore(rate);

           for j = 1:10 % 10x10, 10 bins, each with 10% of the percentile distro
                prctiles_low = prctile(ripples.peakNormedPower,max([j*10-5, 0]));
                prctiles_hi = prctile(ripples.peakNormedPower,min([j*10+5,100]));
                idx = find(ripples.peakNormedPower>prctiles_low & ripples.peakNormedPower<prctiles_hi);
                
                temp = nan(length(idx),1001);
                for id = 1:length(idx) % get the rates for each ripple in this prctile, there are much faster ways to do this...
                    [a b] = min(abs(ripples.peaks(idx(id))-rateTS));
                    if b > 500 & b < length(rate)-500 % keep it in the recording..
                        temp(id,:) = populationZscoredRate(b-500:b+500);
                    end
                end
                prcTileRates_temp(j,:)=nanmean(temp); % average Z-scored PETH around ripples 
                clear temp
           end
           prcTileRates(i,:,:) = prcTileRates_temp;
           imagesc(squeeze(nanmean(prcTileRates))) % show me while you work, computer
           i
           pause(.1)
       end
       end
    end
   clear ripples
   cd ..
end

%% some more plotting
subplot(3,2,1)
imagesc(squeeze(nanmean(prcTileRates)))
ylabel('ripple power prctile')
xlabel('time')

subplot(3,2,2)
plot(nanmean(nanmean(prcTileRates(:,:,450:550),3))) % +/- 50 milliseconds
hold on
plot(nanmean(nanmean(prcTileRates(:,:,400:600),3))) % +/- 100 milliseconds


subplot(3,2,3)
hold on