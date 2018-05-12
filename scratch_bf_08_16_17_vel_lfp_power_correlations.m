load('bf_dataset.mat', 'file')
load('S_markers_times.mat')
load('bf_dataset.mat', 'positions')
load('lfp_trial_avgs.mat', 'lfp_all_recs')
load('lfp_trial_avgs.mat', 'list')
power = dir('LFP_ratemaps_both/*mat');
mark = zeros(length(file),1);

for i=1:length(file)
for j=1:length(power)
if strcmp(file(i).name(1:3),power(j).name(1:3)) & strcmp(file(i).name(8:9),power(j).name(9:10)) 
for l = 1:36
if strcmp(file(i).name(1:3),list(l).name(1:3)) & strcmp(file(i).name(8:9),list(l).name(9:10))
   mark(i) = 1;
end
end
if mark(i) == 1
power_data = load(['LFP_ratemaps_both/' power(j).name]);
for t=1:size(power_data.temp_wave,2)
% [blah start] = min(abs(positions{i}(:,1)-S_markers_times{i}(t,1)));
% [blah stop] = min(abs(positions{i}(:,1)-S_markers_times{i}(t,end)));
% pos = positions{i}(start:stop,2:5);
% pos(pos==0)=nan;
% vel = makeLength(fillmissing(sqrt(sum(diff(pos').^2)),'linear'),540);
for f = 1:150
    
% [a co{j}(t,f)] = corr(vel',squeeze(power_data.temp_wave(1,t,f,:)));

sk{j}(t,f) = skewness(squeeze(power_data.temp_wave(1,t,f,:)));
kurt{j}(t,f) = kurtosis(squeeze(power_data.temp_wave(1,t,f,:)));
st{j}(t,f) = std(squeeze(power_data.temp_wave(1,t,f,:)));
[dip{j}(t,f),pval]=hartigansdipsigniftest(squeeze(power_data.temp_wave(1,t,f,:)),100);

% [a co{j}(t,f,:)] = crosscorr(squeeze(power_data.temp_wave(1,t,f,:)),squeeze(power_data.temp_wave(1,t,f,:)),200);
% for iter = 1:100
%     rand_trial = ceil(rand*size(power_data.temp_wave,2));
% [a co_shuffle{j}(t,f,iter)] = corr(vel',squeeze(power_data.temp_wave(1,rand_trial,f,:)));
% end
end
end
ku(j,:) = median(kurt{j});
end
i
end
end
end


