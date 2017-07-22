% this code attempts to compile phase/rate coding data with peer prediciotn assemblies t
ass = [];
noAss = [];
stren = [];
devs = [];
noAss_phaseonly = [];
ass_phaseonly = [];
d  = dir('*201*');

for i=[1:length(d)]
cd(d(i).name) 
spikes = bz_GetSpikes;
b = dir('*.behavior.mat');
load(b.name);
nBins = length(behavior.events.map{1}.x);
xml = LoadParameters;
if exist('assembliesCrossRegionData_w_theta_sin_cos_coord_vel.mat') & exist([xml.FileName '.positionDecodingGLM_box.cellinfo.mat']) ...

    load([xml.FileName '.positionDecodingGLM_gaussian.cellinfo.mat'])
    load('assembliesCrossRegionData_w_theta_sin_cos_coord_vel.mat');
    positionDecodingGLM=positionDecodingGLM_gaussian
    if isfield(positionDecodingGLM,'dateRun') & length(pairs)>1 & exist('dev')==1
    conditions = length(unique(positionDecodingGLM.results{1}.condition));
    for cell =1:length(positionDecodingGLM.results)
%         t_rate = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_rate',...
%             'GroupingVariables',{'tau','condition'});
%         t_phase = varfun(@mean,positionDecodingGLM.results{cell},'InputVariables','mse_phase_all',...
%             'GroupingVariables',{'tau','condition'});
%         tab = join(t_rate,t_phase);
        for cond = 1:conditions
            if cond <= length(dev)
            % grab phase/rate coding variables
            rows = find(positionDecodingGLM.results{cell}.condition==cond);
            first500ms = find(ismember(positionDecodingGLM.results{cell}.tau(rows),1:4000));
            
            min_mse_rate(cell,cond) = sqrt(min(positionDecodingGLM.results{cell}.mse_rate(rows(first500ms))))./nBins;
            min_mse_phase_all(cell,cond) = sqrt(min(positionDecodingGLM.results{cell}.mse_phase_all(rows(first500ms))))./nBins;
            
            max_mse_rate(cell,cond) = sqrt(max(positionDecodingGLM.results{cell}.mse_rate(rows(first500ms))))./nBins;
            max_mse_phase_all(cell,cond) = sqrt(max(positionDecodingGLM.results{cell}.mse_phase_all(rows(first500ms))))./nBins;
            end
        end
    end
  
        for cond = 1:conditions
            if cond <= length(dev)
            % now grab assemblies
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
            
            if imp > 7 & b > 7 & b < 75 & zerolag < 1.2 & mean(dev{cond}(:,pair))>150
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
             f=[];
            for t=1:length(spikes.times)
                if strcmp(spikes.region{t},'ls')
                    f = [f;t];
                end
            end
            if ~isempty(p)
                ass = [ass; (min_mse_phase_all(p(:,1),cond)-min_mse_rate(p(:,1),cond)) ./ (min_mse_phase_all(p(:,1),cond)+min_mse_rate(p(:,1),cond))];
                ass_phaseonly = [ass_phaseonly; (min_mse_phase_all(p(:,1),cond)) ];
                ff = f(~ismember(f,p(:,1)));
            else
                ff = f;
            end
            noAss = [noAss; (min_mse_phase_all(ff,cond)-min_mse_rate(ff,cond)) ./ (min_mse_phase_all(ff,cond)+min_mse_rate(ff,cond))];        
            noAss_phaseonly = [noAss_phaseonly; (min_mse_phase_all(ff,cond)) ];        
            
            end
            end

    end
        clear pairs dev devControl
if ~isempty(stren)
subplot(2,2,1)
histogram(ass,[-1:.02:1],'Normalization','pdf','FaceColor','g')
hold on
histogram(noAss,[-1:.02:1],'Normalization','pdf','FaceColor','r')


xlabel('phase coding only             neither/both             rate code only')
set(gca,'yscale','log')
ylabel('probability')
title('phase coding strength for gamma assemblies or no assembly')
hold off
subplot(2,2,2)
scatter(stren(:,2),stren(:,5),'.')
xlabel('improvement val')
ylabel('phase-rate comparison')

subplot(2,2,3);
histogram([ass_phaseonly; noAss_phaseonly],[50],'Normalization','pdf','FaceColor','k')
hold on
histogram(ass_phaseonly,50,'Normalization','pdf','FaceColor','g')
histogram(noAss_phaseonly,50,'Normalization','pdf','FaceColor','r')
xlabel('phase coding                               no phase coding')
set(gca,'yscale','log')
ylabel('probability')
hold off

subplot(2,2,4)
scatter(stren(:,1),stren(:,5),'.')
xlabel('improvement val')
ylabel('phase coding val')
title(d(i).name)
pause(.1)
end
end
% 
cd /home/david/datasets/lsDataset/
end




























