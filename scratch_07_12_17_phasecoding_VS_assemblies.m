% this code attempts to compile phase/rate coding data with peer prediciotn assemblies t
ass = [];
noAss = [];
stren = [];
devs = [];
noAss_phaseonly = [];
ass_phaseonly = [];
phaseonly_chance = [];
ls_rec =[]; hpc_rec = [];
ass_rate=[];
noAss_rate =[];
d  = dir('*201*');

for i=1:length(d)
    figure(1)
    clf
%     
% ass = [];
% noAss = [];
% stren = [];
% devs = [];
% noAss_phaseonly = [];
% ass_phaseonly = [];

cd(d(i).name) 
spikes = bz_GetSpikes;
load([d(i).name '.behavior.mat'],'behavior')
b = dir('*.behavior.mat');
load(b.name);
nBins = length(behavior.events.map{1}.x);
xml = LoadParameters;
if exist('assembliesCrossRegion_split_w_theta.mat') & exist([xml.FileName '.positionDecodingGLM_binnedspace_box.cellinfo.mat']) ...

    load([xml.FileName '.positionDecodingGLM_binnedspace_box.cellinfo.mat'])
    load('assembliesCrossRegion_split_w_theta.mat','dev*','pairs');
    positionDecodingGLM=positionDecodingGLM_binnedspace_box;
    if isfield(positionDecodingGLM,'dateRun') & length(pairs)>1 & exist('dev')==1
    conditions = length(unique(positionDecodingGLM.results{1}.condition));
    for cell =1:length(positionDecodingGLM.results)
        t_rate = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_rate',...
            'GroupingVariables',{'tau','condition'});
        t_phase = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_phase_all',...
            'GroupingVariables',{'tau','condition'});
        t_chance = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_chance',...
            'GroupingVariables',{'tau','condition'});
        tab = join(join(t_rate,t_phase),t_chance);
        for cond = 1:conditions
            if cond <= length(dev)  && sum(behavior.events.trialConditions==cond) > 7
            % grab phase/rate coding variables
            rows = find(tab.condition==cond);
            first500ms = find(ismember(tab.tau(rows),1:nBins/4));
            
%             min_mse_rate(cell,cond) = sqrt(min(tab.mean_mse_rate(rows(first500ms))))./nBins;
%             min_mse_phase_all(cell,cond) = sqrt(min(tab.mean_mse_phase_all(rows(first500ms))))./nBins;
            min_mse_rate(cell,cond) = (min(tab.mean_mse_rate(rows(first500ms)))./mean(tab.mean_mse_chance(rows(first500ms))));
            min_mse_phase_all(cell,cond) = (min(tab.mean_mse_phase_all(rows(first500ms)))./mean(tab.mean_mse_chance(rows(first500ms))));
            
            
            min_mse_chance(cell,cond) = (min(tab.mean_mse_chance(rows(first500ms))))./mean(tab.mean_mse_chance(rows(first500ms)));
            
            max_mse_rate(cell,cond) = sqrt(max(tab.mean_mse_rate(rows(first500ms))))./nBins;
            max_mse_phase_all(cell,cond) = sqrt(max(tab.mean_mse_phase_all(rows(first500ms))))./nBins;
            end
        end
    end
  
        pairCount = 0;
        h=[];
        z=[];
        ii=[];
        p=[];
        for cond = 1:conditions
            if cond <= length(dev) && sum(behavior.events.trialConditions==cond) > 7 % check that assemblies have run and there are enough trials
            % now grab assemblies
            for pair = 1:size(dev{cond},2)
            [a b] =  min(dev{cond}(:,pair));
            [aa bb] = min(mean(devControl{cond}(:,pair,:),3));
            imp = (a-mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(:,pair,:),3)));
            imp2 = a ./ max(mean(mean(devControl{cond}(:,pair,:),3)));
            zerolag = (min(dev{cond}(1:6,pair)) - mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(1,pair,:),3)));
            if zerolag < 1 
                zerolag = 1;
            end
            if imp > 4.5 & b > 7 & b < 150 &  zerolag < 1.2 & mean(dev{cond}(:,pair))>100
                p = [p; pairs(pair,:)];
                h = [h; imp];
                ii = [ii;b];
                z=[z;zerolag];
            end
            pairCount = 1 + pairCount;
            end
            end
        end
        for cond = 1:conditions
            if cond <= length(dev) && sum(behavior.events.trialConditions==cond) > 10 % check that assemblies have run and there are enough trials
            % now grab assemblies
%             if length(h) > 8
                pairCount = 0;
                h=[];
                z=[];
                ii=[];
                p=[];
                for pair = 1:size(dev{cond},2)
                [a b] =  min(dev{cond}(:,pair));
                [aa bb] = min(mean(devControl{cond}(:,pair,:),3));
                imp = (a-mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(:,pair,:),3)));
                imp2 = a ./ max(mean(mean(devControl{cond}(:,pair,:),3)));
                zerolag = (min(dev{cond}(1:6,pair)) - mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(1,pair,:),3)));
                if zerolag < 1 
                zerolag = 1;
                end

                if imp > 4.5 & b > 7 & b < 150 &  zerolag < 1.2 & mean(dev{cond}(:,pair))>100
                p = [p; pairs(pair,:)];
                h = [h; imp];
                ii = [ii;b];
                z=[z;zerolag];
                stren = [stren; imp imp2 b (min_mse_phase_all(pairs(pair,1),cond)-min_mse_rate(pairs(pair,1),cond)) ...
                ./ (min_mse_phase_all(pairs(pair,1),cond)+min_mse_rate(pairs(pair,1),cond)) ...
                (min_mse_phase_all(pairs(pair,1),cond))];
                devs = [devs;dev{cond}(:,pair)'];
                end
                pairCount = 1 + pairCount;
                end 
%             end
             f=[];
            for t=1:length(spikes.times)
                if strcmp(spikes.region{t},'ls')
                    f = [f;t];
                end
            end
            if ~isempty(p)
                ass = [ass; (min_mse_phase_all(p(:,1),cond)-min_mse_rate(p(:,1),cond)) ./ (min_mse_phase_all(p(:,1),cond)+min_mse_rate(p(:,1),cond))];
                ass_phaseonly = [ass_phaseonly; (min_mse_phase_all(p(:,1),cond)) ];
                ass_rate = [ass_rate; (min_mse_rate(p(:,1),cond)) ];
                phaseonly_chance = [phaseonly_chance; min_mse_chance(p(:,1),cond)];
                ff = f(~ismember(f,p(:,1)));
                ls_rec = [ls_rec;  repmat(i,size(p,1),1)];
            else
                ff = f;
            end
            noAss = [noAss; (min_mse_phase_all(ff,cond)-min_mse_rate(ff,cond)) ./ (min_mse_phase_all(ff,cond)+min_mse_rate(ff,cond))];        
            noAss_phaseonly = [noAss_phaseonly; (min_mse_phase_all(ff,cond)) ];  
            noAss_rate = [noAss_rate; (min_mse_rate(ff,cond)) ];  
            phaseonly_chance = [phaseonly_chance; min_mse_chance(ff,cond)];
            ls_rec = [ls_rec; repmat(i,length(ff),1)];
            end
            end

    end
       
if ~isempty(stren)
    noAss_phaseonly(noAss_phaseonly>1)=1;
    noAss(noAss>1)=1;
    ass_phaseonly(ass_phaseonly>1)=1;
    ass(ass>1)=1;
    
subplot(2,2,1)
histogram(ass,[-1:.02:1],'Normalization','pdf','FaceColor','g')
hold on
histogram(noAss,[-1:.02:1],'Normalization','pdf','FaceColor','r')
line([mean(ass) mean(ass)],[0 40],'color','g')
line([mean(noAss) mean(noAss)],[0 40],'color','r')
xlabel('phase coding only             neither/both             rate code only')
% set(gca,'yscale','log')
ylabel('probability')
title(length(pairs))
% title('phase coding strength for gamma assemblies or no assembly')
hold off
subplot(2,2,2)
scatter(stren(:,2),stren(:,5),'.')
xlabel('improvement val')
ylabel('phase-rate comparison')

subplot(2,2,3);
% histogram([ass_phaseonly; noAss_phaseonly],[50],'Normalization','pdf','FaceColor','k')
hold on
histogram(ass_phaseonly,0:.01:1,'Normalization','pdf','FaceColor','g')
histogram(noAss_phaseonly,0:.01:1,'Normalization','pdf','FaceColor','r')
xlabel('phase coding                               no phase coding')
% set(gca,'yscale','log')
line([mean(ass_phaseonly) mean(ass_phaseonly)],[0 25],'color','g')
line([mean(noAss_phaseonly) mean(noAss_phaseonly)],[0 25],'color','r')
ylabel('probability')
title([num2str(size(stren,1)) ' assemblies detected'])
hold off

subplot(2,2,4)
scatter(stren(:,1),stren(:,5),'.')
xlabel('improvement val')
ylabel('phase coding val')
title(d(i).name)
[a b] = kstest2(ass_phaseonly,noAss_phaseonly);
subplot(2,2,2)
title(b);


% figure(100)
% subplot(2,2,1)
% scatter((size(stren,1)),b,'.k')
% set(gca,'yscale','log')
% title('# assemblies vs pval')
% hold on
% subplot(2,2,2)
% scatter(mean(ass_phaseonly)-mean(noAss_phaseonly),b,'.k')
% title('phase coding diff vs pval')
% set(gca,'yscale','log')
% hold on
% subplot(2,2,3)
% scatter(mean(ass_phaseonly)-mean(noAss_phaseonly),(size(stren,1)),'.k')
% title('phase coding diff vs # of assemblies')
% hold on

pause(.01)
end
 clear pairs dev devControl
end
% 

    cd('/home/david/datasets/lsDataset')
% cd('/home/david/datasets/lsDataset')
end










figure
clear a n c ap np cp 
for iter = 1:1000
r = randperm(length(ass_phaseonly));
a(iter,:) = ass_phaseonly(r(1:round(prctile(1:length(ass_phaseonly),60))));
r = randperm(length(noAss_phaseonly));
n(iter,:) = noAss_phaseonly(r(1:round(prctile(1:length(noAss_phaseonly),60))));
r = randperm(length(phaseonly_chance));
c(iter,:) = phaseonly_chance(r(1:round(prctile(1:length(phaseonly_chance),60))));
for j = 1:100
ap(iter,j) = prctile(a(iter,:),j);
np(iter,j) = prctile(n(iter,:),j);
cp(iter,j) = prctile(c(iter,:),j);
end
boundedline(1:100,mean(ap),3*std(ap),'g');
boundedline(1:100,mean(np),3*std(np),'r');
boundedline(1:100,mean(cp),3*std(cp),'k');
axis([0 100 .5 1])
title(iter)
pause(.0001);
end




