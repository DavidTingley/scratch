d = dir('*/*firingMaps*');
b = dir('*/*behavior*');
l = 1; h=1;

for i=1:length(d)
    load([d(i).folder '/' d(i).name]);
    cd(d(i).folder)
    b = dir('*behavior*');
    load([b(1).name]);
    cd /home/david/datasets/lsDataset
    if strcmp(behavior.description,'wheel alternation')
    [binnedfiringMaps.phaseMaps] = bz_phaseMap2Bins(firingMaps.phaseMaps,firingMaps.rateMaps,behavior);
    for cell = 1:length(firingMaps.UID)
        for condition = 1:length(firingMaps.rateMaps)
           if strcmp(firingMaps.region{cell},'ls') & size(firingMaps.rateMaps{condition},2) > 10
            l_maps(l,:) = makeLength(squeeze(mean(firingMaps.rateMaps{condition}(cell,:,:),2)),200);
            l_phase_maps(l,:) = makeLength(squeeze(circ_mean(squeeze(binnedfiringMaps.phaseMaps{condition}(cell,:,:)))),200);
            l=l+1;
           elseif size(firingMaps.rateMaps{condition},2) > 10
            h_maps(h,:) = makeLength(squeeze(mean(firingMaps.rateMaps{condition}(cell,:,:),2)),200);
            h_phase_maps(h,:) = makeLength(squeeze(circ_mean(squeeze(binnedfiringMaps.phaseMaps{condition}(cell,:,:)))),200);
            h=h+1;
           end
        end
    end  
    end
end

for i=1:length(h_phase_maps)
h_phase_maps_smooth(i,:) = circ_smoothTS(h_phase_maps(i,:),30,'method','median','exclude',0);
end
for i=1:length(l_phase_maps)
l_phase_maps_smooth(i,:) = circ_smoothTS(l_phase_maps(i,:),30,'method','median','exclude',0);
end 
hp = circ_mean(h_phase_maps_smooth);
lp = circ_mean(l_phase_maps_smooth);



plot(hp)
figure
plot(lp,'g')