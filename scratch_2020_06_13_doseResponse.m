%% turn this into a dose response comparison for ripple rate AND theta
th = 30;


cgm_files = dir('cgm*mat'); % opto animals

inc = [1 3 4 5 6 8];
keep{1} = 1:1400;
keep{2} = [];
keep{3} = 945:2500;
keep{4} = 1:2000;
keep{5} = 1:1700;
keep{6} = 1190:3120;
keep{7} = [];
keep{8} = 1:1500;


for ff = 1:8 % animals
load([cgm_files(ff).folder '/' cgm_files(ff).name],'stimRate','count','theta_*','states','spSlope','theta_resid','emgSig','idx','isa','isig_levels')
c=1;
stimRate(end)=nan;
id = intersect(keep{ff},idx);

if ff ==1
    stimRate = stimRate/270*200;
elseif ff == 3
    stimRate = stimRate/305*200;  
elseif ff == 4
    stimRate = stimRate/240*200;
elseif ff == 6 
    stimRate = stimRate/380*200;
elseif ff == 8 
    stimRate = stimRate/600*200;
end


d = (stimRate);
d(d==0)=nan;
    for i=1:100
        thresh = prctile(d(idx),max([0, i-15]));
        thresh_hi = prctile(d(idx),min([100, i+15]));
        f = find(d>thresh & d<=thresh_hi);
        f = (intersect(f,idx));
        jj=[];
        for j=1:length(f)
        if f(j)>20 & f(j) < length(d)-20
        jj(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
        end
        end
        if ~isempty(jj)
        doseResp(ff,c,:) = nanmean(jj,1);clear jj;
        end
        c=1+c;
    end
    for i=1:400
       f = find(d>i-th & d<i+th);
%        f = find(d==i);
       f = (intersect(f,idx));
       jj=[];
       for j=1:length(f)
        if f(j)>20 & f(j) < length(d)-20
            jj(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
        end
        end
        if ~isempty(jj)
            doseResp_rate(ff,i,:) = nanmean(jj,1);clear jj;
        else
            doseResp_rate(ff,i,:) = nan(1,40);
        end 
    end
    for it = 1:100
        for i=1:200
           offset = randi(100)-50;
           f = find(circshift(d,offset)==i);
           f = (intersect(f,idx));
           jj=[];
           for j=1:length(f)
            if f(j)>20 & f(j) < length(d)-20
                jj(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
            end
            end
            if ~isempty(jj)
                doseResp_rate_shuf(ff,it,i,:) = nanmean(jj,1);clear jj;
            else
                doseResp_rate_shuf(ff,it,i,:) = nan(1,40);
            end
        end
    end
    
end

for i=1:8
s(i,:) = squeeze(mean(doseResp(i,:,23:27),3));
s_rate(i,:) = squeeze(nanmean(doseResp_rate(i,:,23:27),3));
[c(i) p(i)] = corr(s_rate(i,:)',[1:size(s_rate,2)]','rows','complete');
for it = 1:100
   s_rate_shuf(i,it,:)=squeeze(nanmean(doseResp_rate_shuf(i,it,:,23:27),4));
   [cs(i,it),p(i,it)]=corr(squeeze(s_rate_shuf(i,it,:)),[1:size(s_rate_shuf,3)]','rows','complete');
end
end




%%%%%%%%%%%%%%%%%%%%

cd ..
% cgm_files = dir('_*mat'); % no stim animals
orig_files = dir('_*mat');
new_files = dir('revision/*mat');

cgm_files = [orig_files;new_files];

for ff = 1:length(cgm_files) % animals
load([cgm_files(ff).folder '/' cgm_files(ff).name],'stimRate','count','theta_*','states','spSlope','theta_resid','emgSig','idx','isa','isig_levels')
c=1;
d = (count);
d(d==0)=nan;
for i=1:95
    thresh = prctile(d(idx),max([0, i-5]));
    thresh_hi = prctile(d(idx),min([100, i+5]));
    f = find(d>thresh & d<=thresh_hi);
    f = intersect(f,idx);jj=[];
    for j=1:length(f)
    if f(j)>20 & f(j) < length(d)-20
    jj(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
    end
    end
    if ~isempty(jj)
    doseResp(ff,c,:) = nanmean(jj,1);clear jj;
    end
    c=1+c;
end
for i=1:300
   f = find(d>i-th & d<i+th);
%    f = find(d==i);
   f = intersect(f,idx);
   jj=[];
   for j=1:length(f)
    if f(j)>20 & f(j) < length(d)-20
        jj(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
    end
    end
    if ~isempty(jj)
        doseResp_rate(ff,i,:) = nanmean(jj,1);clear jj;
    else
        doseResp_rate(ff,i,:) = nan;
    end 
end

for it = 1:100
    for i=1:200
       offset = randi(100)-50;
       f = find(circshift(d,offset)==i);
       f = intersect(f,idx);
       jj=[];
       for j=1:length(f)
        if f(j)>20 & f(j) < length(d)-20
            jj(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
        end
        end
        if ~isempty(jj)
            doseResp_rate_shuf(ff,it,i,:) = nanmean(jj,1);clear jj;
        else
            doseResp_rate_shuf(ff,it,i,:) = nan(1,40);
        end
    end
end
    
end

for i=1:length(cgm_files)
s(i,:) = squeeze(mean(doseResp(i,:,21:23),3));
s_rate(i,:) = squeeze(nanmean(doseResp_rate(i,:,21:23),3));
[c(i) p(i)] = corr(s_rate(i,:)',[1:size(s_rate,2)]','rows','complete');
for it = 1:100
   s_rate_shuf(i,it,:)=squeeze(nanmean(doseResp_rate_shuf(i,it,:,21:23),4));
   [cs(i,it),p(i,it)]=corr(squeeze(s_rate_shuf(i,it,:)),[1:size(s_rate_shuf,3)]','rows','complete');
end
end
s_rate(s_rate==0)=nan;
s(s==0)=nan;
% 
