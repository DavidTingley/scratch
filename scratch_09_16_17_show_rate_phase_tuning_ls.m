for i=76:113
for c = 1:9
figure(1); subplot(2,5,c);figure(2);subplot(2,5,c)
for spk = 1:size(firingMaps.phaseMaps{c}{i},1)
if spk < 300
trial = firingMaps.phaseMaps{c}{i}(spk,2);
bin = firingMaps.phaseMaps{c}{i}(spk,1);
figure(1); scatter(firingMaps.phaseMaps{c}{i}(spk,end),firingMaps.rateMaps{c}(44,trial,bin),'.k')
hold on
figure(2)
scatter(firingMaps.rateMaps{c}(i,trial,bin),firingMaps.rateMaps{c}(44,trial,bin),'.r')
hold on
end
end
end
pause;figure(1); clf; figure(2); clf
end