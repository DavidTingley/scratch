dd = dir('*201*')

for i=1:length(dd)
    if dd(i).isdir
   cd(dd(i).name)
   if exist('behav.mat') & isempty(dir('*behavior.mat'))
       load('behav.mat','trials','pos','map*')
       if size(pos,2) == 10
          Process_ConvertOptitrack2Pos(dd(i).name)
          b = dir('*pos');
          pos = importdata(b.name);
            for i=1:length(trials)
            for j=1:length(trials{i})
            [a start] = min(abs(pos(:,1)-trials{i}{j}(1,1)));
            [a stop] = min(abs(pos(:,1)-trials{i}{j}(end,1)));
            trials_new{i}{j} = pos(start:stop,:);
            end
            end
            trials = trials_new;clear trials_new;clear map mapping
            c=1;bins = 200;
            for tt = 1:length(trials)
%             if d == 5  %% only for rec 21st
%                 bins = 100;
%             elseif d == 4
%                 bins = 100;
%             else
%                 bins = 80;
%             end
    
    hold off
    map{tt}=[];
    t_conc=[];
    for t = 1:length(trials{tt})
        t_conc = [trials{tt}{t}(:,:),20*(trials{tt}{t}(:,1)-trials{tt}{t}(1,1))];
    if length(t_conc)>=bins
        while length(t_conc)>bins+1
        di = pdist(t_conc);
        s = squareform(di);
        s(find(eye(size(s))))=nan;
        [a b] = min(s(:));
        [coords blah] = find(s==a);
        t_conc(coords(1),:) = (t_conc(coords(1),:)+t_conc(coords(2),:))./2;
        t_conc(coords(2),:) = [];
        % debug
%         scatter(t_conc(:,1),t_conc(:,2))
%         pause(.01)
        end
    t_conc_all(t,:,:) = t_conc;
    end
    end
    if length(trials{tt})>0
    map{tt} = squeeze(median(t_conc_all(:,:,:),1));

    if ~isempty(trials{tt})
    subplot(8,4,c)
    c=1+c
    for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
    scatter(trials{tt}{t}(:,1),trials{tt}{t}(:,2),'.k')
    hold on
    axis([0 550 0 550])
    end
    for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
    scatter(trials{tt}{t}(1,1),trials{tt}{t}(1,2),'.g')
    scatter(trials{tt}{t}(end,1),trials{tt}{t}(end,2),'.r')
    hold on
    
    axis([0 550 0 550])
    end
    end
    end
    clear t_conc_all
%     for t = 1:length(trials{tt})
%         scatter((trials{tt}{t}(:,1)),(trials{tt}{t}(:,2)),'.r')
%         hold on
%         en
%         scatter((trials{tt}{t}(:,3)),(trials{tt}{t}(:,4)),'.b')
%         scatter((trials{tt}{t}(1,3)),(trials{tt}{t}(1,4)),'.g')
%     end
%     if ~isempty(trials{tt})
%     g = ginput();
%     ggx=[];ggy=[];
%     for t =1:length(g)-1whps
%         dd = abs(g(t,1)-g(t+1,1)) + abs(g(t,2)-g(t+1,2));
%         ggx = [ggx,linspace(g(t,1),g(t+1,1),round(dd))];
%     
%         ggy = [ggy,linspace(g(t,2),g(t+1,2),round(dd))];
%     end
%     g = [ggx;ggy]';
%     
%     step = (max(g(:,2))-min(g(:,2)))./bins;
%     step2 = (max(g(:,2))-min(g(:,2)))./length(g(:,2));
%     ww = smooth(interp1(min(g(:,2)):step2+.0001:max(g(:,2)),g(:,2),min(g(:,2)):step:max(g(:,2))+.01,'linear','extrap'),3);
%     step = (max(g(:,1))-min(g(:,1)))./bins;
%     step2 = (max(g(:,1))-min(g(:,1)))./length(g(:,1));
%     w = smooth(interp1(min(g(:,1)):step2+.0001:max(g(:,1)),g(:,1),min(g(:,1)):step:max(g(:,1))+.01,'linear','extrap'),3);
%     plot(w,ww,'k')
%     disp('finding mapping...')
%     map{tt} = [w,ww];
    for t =1:length(trials{tt})  % all trial types (rotations)
        for p = 1:length(trials{tt}{t})
            [a b] = min(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
                trials{tt}{t}(p,8)-map{tt}(:,8),...
                trials{tt}{t}(p,9)-map{tt}(:,9),...
                trials{tt}{t}(p,10)-map{tt}(:,10),...
                (trials{tt}{t}(p,1)-trials{tt}{t}(1,1))*50-map{tt}(:,1),...  % penalty for time differences
                40*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'));     % penalty for order differences
            mapping{tt}{t}(p,:) = [map{tt}(b,1:end) b trials{tt}{t}(p,1)];
%             plot(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
%                 trials{tt}{t}(p,2)-map{tt}(:,2),...
%                 trials{tt}{t}(p,3)-map{tt}(:,3),...
%                 trials{tt}{t}(p,4)-map{tt}(:,4),...
%                 (trials{tt}{t}(p,5)-trials{tt}{t}(1,5))*20-map{tt}(:,5),...  % penalty for time differences
%                 50*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'))
%             pause
        end
    end
            end % regenerate mapping data
       end
       behavior = pos2behav(pos,'optitrack','trials',trials,'mapping',mapping,'map',map,'behavType','wheel alternation');
       save([dd(i).name '.behavior.mat'],'behavior');
   end
   if exist('behav_alt.mat') & isempty(dir('*behavior.mat'))
       load('behav_alt.mat','trials','pos','map*')
       if size(pos,2) == 10
          Process_ConvertOptitrack2Pos(dd(i).name)
          b = dir('*pos');
          pos = importdata(b.name);
            for i=1:length(trials)
            for j=1:length(trials{i})
            [a start] = min(abs(pos(:,1)-trials{i}{j}(1,1)));
            [a stop] = min(abs(pos(:,1)-trials{i}{j}(end,1)));
            trials_new{i}{j} = pos(start:stop,:);
            end
            end
            trials = trials_new;clear trials_new;clear map mapping
            c=1;bins = 200;
            for tt = 1:length(trials)
%             if d == 5  %% only for rec 21st
%                 bins = 100;
%             elseif d == 4
%                 bins = 100;
%             else
%                 bins = 80;
%             end
    
    hold off
    map{tt}=[];
    t_conc=[];
    for t = 1:length(trials{tt})
        t_conc = [trials{tt}{t}(:,:),20*(trials{tt}{t}(:,1)-trials{tt}{t}(1,1))];
    if length(t_conc)>=bins
        while length(t_conc)>bins+1
        di = pdist(t_conc);
        s = squareform(di);
        s(find(eye(size(s))))=nan;
        [a b] = min(s(:));
        [coords blah] = find(s==a);
        t_conc(coords(1),:) = (t_conc(coords(1),:)+t_conc(coords(2),:))./2;
        t_conc(coords(2),:) = [];
        % debug
%         scatter(t_conc(:,1),t_conc(:,2))
%         pause(.01)
        end
    t_conc_all(t,:,:) = t_conc;
    end
    end
    if length(trials{tt})>0
    map{tt} = squeeze(median(t_conc_all(:,:,:),1));

    if ~isempty(trials{tt})
    subplot(8,4,c)
    c=1+c
    for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
    scatter(trials{tt}{t}(:,1),trials{tt}{t}(:,2),'.k')
    hold on
    axis([0 550 0 550])
    end
    for t = 1:length(trials{tt})
%     scatter(map{tt}(:,1),map{tt}(:,2),'.')
    scatter(trials{tt}{t}(1,1),trials{tt}{t}(1,2),'.g')
    scatter(trials{tt}{t}(end,1),trials{tt}{t}(end,2),'.r')
    hold on
    
    axis([0 550 0 550])
    end
    end
    end
    clear t_conc_all
%     for t = 1:length(trials{tt})
%         scatter((trials{tt}{t}(:,1)),(trials{tt}{t}(:,2)),'.r')
%         hold on
%         en
%         scatter((trials{tt}{t}(:,3)),(trials{tt}{t}(:,4)),'.b')
%         scatter((trials{tt}{t}(1,3)),(trials{tt}{t}(1,4)),'.g')
%     end
%     if ~isempty(trials{tt})
%     g = ginput();
%     ggx=[];ggy=[];
%     for t =1:length(g)-1whps
%         dd = abs(g(t,1)-g(t+1,1)) + abs(g(t,2)-g(t+1,2));
%         ggx = [ggx,linspace(g(t,1),g(t+1,1),round(dd))];
%     
%         ggy = [ggy,linspace(g(t,2),g(t+1,2),round(dd))];
%     end
%     g = [ggx;ggy]';
%     
%     step = (max(g(:,2))-min(g(:,2)))./bins;
%     step2 = (max(g(:,2))-min(g(:,2)))./length(g(:,2));
%     ww = smooth(interp1(min(g(:,2)):step2+.0001:max(g(:,2)),g(:,2),min(g(:,2)):step:max(g(:,2))+.01,'linear','extrap'),3);
%     step = (max(g(:,1))-min(g(:,1)))./bins;
%     step2 = (max(g(:,1))-min(g(:,1)))./length(g(:,1));
%     w = smooth(interp1(min(g(:,1)):step2+.0001:max(g(:,1)),g(:,1),min(g(:,1)):step:max(g(:,1))+.01,'linear','extrap'),3);
%     plot(w,ww,'k')
%     disp('finding mapping...')
%     map{tt} = [w,ww];
    for t =1:length(trials{tt})  % all trial types (rotations)
        for p = 1:length(trials{tt}{t})
            [a b] = min(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
                trials{tt}{t}(p,8)-map{tt}(:,8),...
                trials{tt}{t}(p,9)-map{tt}(:,9),...
                trials{tt}{t}(p,10)-map{tt}(:,10),...
                (trials{tt}{t}(p,1)-trials{tt}{t}(1,1))*50-map{tt}(:,1),...  % penalty for time differences
                40*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'));     % penalty for order differences
            mapping{tt}{t}(p,:) = [map{tt}(b,1:end) b trials{tt}{t}(p,1)];
%             plot(nansum(abs([trials{tt}{t}(p,1)-map{tt}(:,1),...
%                 trials{tt}{t}(p,2)-map{tt}(:,2),...
%                 trials{tt}{t}(p,3)-map{tt}(:,3),...
%                 trials{tt}{t}(p,4)-map{tt}(:,4),...
%                 (trials{tt}{t}(p,5)-trials{tt}{t}(1,5))*20-map{tt}(:,5),...  % penalty for time differences
%                 50*(p./length(trials{tt}{t})*length(map{tt}) - (1:length(map{tt})))'])'))
%             pause
        end
    end
            end % regenerate mapping data
       end
       behavior = pos2behav(pos,'optitrack','trials',trials,'mapping',mapping,'map',map,'behavType','wheel alternation');
       save([dd(i).name '.behavior.mat'],'behavior');
   end
   cd ..   
    end
end