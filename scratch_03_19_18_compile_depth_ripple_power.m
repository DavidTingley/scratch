clear all
recCount = 1;
% cd /home/david/datasets/ripples_LS
cd E:\datasets\ripples_LS
recordingList  = dir('*201*');
animals = {'DT1','DT2','DT5','DT7','DT8','DT9'};
colors = {'c','r','m',[.91 .41 .17],'b','g'};
offsets = [-1200 -2800 -800 0 1000 200];

for i=1:length(recordingList)
   cd(recordingList(i).name) 
   sessionInfo = bz_getSessionInfo;
   if exist([sessionInfo.FileName '.LSRipples.events.mat'])
%    if strcmp(sessionInfo.animal,'DT7')
       load([sessionInfo.FileName '.LSRipples.events.mat'])
       states = bz_LoadStates(pwd,'SleepState');
       nremTime = sum([1; diff(states.ints.NREMstate')']);
       lfp = bz_GetLFP(1);
   if length(ripples.peaks) > 15
%        avgWaveSpec{i} = zeros(length(ripples.peaks),100,625,'single');
       avgRippleBand{i} = zeros(length(ripples.peaks),625,'single');
       ripCount{i} = 0;
       
       for ripple = 1:length(ripples.peaks)
          lfp = bz_GetLFP(ripples.rippleChan,'restrict',...
              [ripples.peaks(ripple)-.25001 ripples.peaks(ripple)+.25001]);
%            for ch = 1%1:size(lfp.data,2)
               wavespec = bz_WaveSpec(lfp.data,'frange',[1 200],'ncyc',3,...
                   'nfreqs',100,'space','lin','samplingRate',1250);
%                avgWaveSpec{i}(ripple,:,:) = abs(wavespec.data');
               avgRippleBand{i}(ripple,:) = (mean(abs(wavespec.data(:,60:100)')));
               ripCount{i} = ripCount{i} + 1;
%            end
       end      
%        avgWaveSpec{i} =  squeeze(median(avgWaveSpec{i},1));
       avgRippleBand{i} =  squeeze(median(avgRippleBand{i},1));
       avgRippleBand_z{i} =  zscore(squeeze(median(avgRippleBand{i},1)));
%        nremRips = Restrict(ripples.timestamps,double(states.ints.NREMstate));
       ripRate(i) = length(ripples.timestamps)./lfp.timestamps(end);
%        NREMTime(i) = nremTime;
       depth(i) = sessionInfo.depth;
       animal{i} = sessionInfo.animal;    
       chan(i) = ripples.rippleChan;
       for ii = 1:6
           ind(ii) = strcmp(animal{i},animals{ii});
       end
%        subplot(3,2,find(ind));
%        plot(max(avgRippleBand{i}(300:325)),-depth(i),'.k')
%        hold on
       
%        pause(.001)
       
%    subplot(4,4,recCount)
%    imagesc(-.5:1/1250:.4999,sessionInfo.spikeGroups.groups{7},flipud(avgRippleBand{i}));
%    title(['day: ' num2str(recCount) ', depth: ' num2str(sessionInfo.depth) ', ripples: ' num2str(ripple)])
   recCount = 1+recCount
   end
   end
%    end
%    pause(.1)
%    cd /home/david/datasets/ripples_LS
cd E:\datasets\ripples_LS
end

for i=1:6
    subplot(2,6,i+6)
    for j=1:length(animal)
        ind(j) = strcmp(animals{i},animal{j});
        if ~isempty(avgRippleBand{j})
        m(j) = max(avgRippleBand{j}(300:325));
        if ind(j) == 1
        plot(max(avgRippleBand{j}(300:325)),-depth(j)-offsets(i),'.k')
        hold on
        end
        end
    end
    idx = find(ind);
    pow = m(idx);
    dep = depth(idx);
    
    [p s]= polyfit(dep,pow,1);
    f{1} = polyval(p,dep);
%     plot(f{1},-dep,'color',colors{i})
    pdf = normpdf(pow,f{1});
    pdf(pdf==0)=nan;
    logL1 = nansum(log(pdf));
    
    
    [p ss] = polyfit(dep,pow,2);
    f{2} = polyval(p,dep);
%     plot(f,-dep,'color',colors{i})
    pdf = normpdf(pow,f{2});
    pdf(pdf==0)=nan;
    logL2 = nansum(log(pdf));
   
    
    [aic,bic] = aicbic([logL1 logL2],[2 1],length(pow));
    [a mod]=min(aic);
    plot(f{mod},-dep-offsets(i),'color',colors{i})
    title([animals{i} ', polyfit: ' num2str(mod)])
    axis([0 1600 -2000 0])
end


% save('/home/david/datasets/ripples_LS/DT7_ripplePower_allChans.mat','-v7.3')