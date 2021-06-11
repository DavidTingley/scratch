cd data
cgm_files = dir('_*.mat');
hourHz = 1/(60*60);
dayHz = 1/(24*60*60);
nyquist = 5./(24*60)/2;

%% plot PSD
for ff=1:length(cgm_files)
    data{ff} = load(cgm_files(ff).name,'isa','absTime','zt','idx','theta*','count','isig_levels','spSlope','states','emgSig','mov','ccc');
end
[b a] = butter(4,[ dayHz ./ 1.5 / nyquist],'high');
for i=1:8
iss = fillmissing([0,diff(data{i}.isig_levels)],'linear');
iss(isnan(data{i}.isig_levels)) = 0;
[pxx(i,:) f] = periodogram(iss(data{i}.idx),[],256,1/(60*5));

iss = fillmissing(data{i}.isig_levels,'linear');
iss = filtfilt(b,a,iss);
iss(isnan(data{i}.isig_levels)) = 0;

[pxx_filt(i,:) f] = periodogram(iss(data{i}.idx),[],256,1/(60*5));

[acg(i,:) lags bounds(i,:)] = autocorr(diff(data{i}.isig_levels(data{i}.idx)),'NumLags',50,'NumSTD',3);
% plot(f,pxx,'k'); p(i,:) = pxx;
% hold onfor
end


cd revision
cgm_files = dir('*.mat');
for ff=1:length(cgm_files)
    data{ff} = load(cgm_files(ff).name,'isa','absTime','zt','idx','theta*','count','isig_levels','spSlope','states','emgSig','mov','ccc');
end
%% plot PSD
[b a] = butter(4,[ dayHz ./ 1.5 / nyquist],'high');
for i=[1:length(cgm_files)]
iss = fillmissing([0,diff(data{i}.isig_levels)],'linear');
iss(isnan(data{i}.isig_levels)) = 0;
if length(data{i}.idx)>50
[pxx(i,:) f] = periodogram(iss(data{i}.idx),[],256,1/(60*5));

iss = fillmissing(data{i}.isig_levels,'linear');
iss = filtfilt(b,a,iss);
iss(isnan(data{i}.isig_levels)) = 0;

[pxx_filt(i,:) f] = periodogram(iss(data{i}.idx),[],256,1/(60*5));

[acg(i+8,:) lags bounds(i+8,:)] = autocorr(diff(data{i}.isig_levels(data{i}.idx)),'NumLags',50,'NumSTD',3);
end
% plot(f,pxx,'k'); p(i,:) = pxx;
% hold onfor
end