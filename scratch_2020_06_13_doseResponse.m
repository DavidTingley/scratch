%% turn this into a dose response comparison for ripple rate AND theta
% cgm_files = dir('cgm*mat'); % opto animals
% 
%   
% for ff = 1:8 % animals
% load(cgm_files(ff).name,'stimRate','count','theta_*','states','spSlope','theta_resid','emgSig','idx','isa','isig_levels')
% c=1;
% d = (stimRate);
% d(d==0)=nan;
% for i=1:95
%     thresh = prctile(d,i);
%     thresh_hi = prctile(d,i+5);
%     f = find(d>thresh & d<=thresh_hi);
%     f = intersect(f,idx);jj=[];
%     for j=1:length(f)
%     if f(j)>20 & f(j) < length(d)-20
%     jj(j,:) = diff(isig_levels(f(j)-20:f(j)+20));
%     end
%     end
%     if ~isempty(jj)
%     doseResp(ff,c,:) = nanmean(jj,1);clear jj;
%     end
%     c=1+c;
%     end
% end
% for i=1:8
% s(i,:) = squeeze(mean(doseResp(i,:,23:27),3));
% [c(i) p(i)] = corr(s(i,:)',[1:size(s,2)]','rows','complete');
% end


cgm_files = dir('_*mat'); % no stim animals

for ff = 1:8 % animals
load(cgm_files(ff).name,'stimRate','count','theta_*','states','spSlope','theta_resid','emgSig','idx','isa','isig_levels')
c=1;
d = (count);
d(d==0)=nan;
for i=1:95
    thresh = prctile(d,i);
    thresh_hi = prctile(d,i+5);
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
for i=1:200
%    f = find(d>i-15 & d<i+15);
   f = find(d==i);
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
            doseResp_rate_shuf(ff,it,i,:) = nan;
        end
end
    
end

for i=1:8
s(i,:) = squeeze(mean(doseResp(i,:,21:23),3));
s_rate(i,:) = squeeze(nanmean(doseResp_rate(i,:,21:23),3));
[c(i) p(i)] = corr(s(i,:)',[1:size(s,2)]','rows','complete');
end
s_rate(s_rate==0)=nan;
s(s==0)=nan;

