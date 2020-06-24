%% turn this into a dose response comparison for ripple rate AND theta
 cgm_files = dir('_*mat');
for ff = 1:8 % animals
load(cgm_files(ff).name,'count','theta_*','states','spSlope','theta_resid','emgSig','idx','isa','isig_levels')
c=1;
d = (emgSig);
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
end
for i=1:8
s(i,:) = squeeze(mean(doseResp(i,:,23:27),3));
[c(i) p(i)] = corr(s(i,:)',[1:size(s,2)]','rows','complete');
end


