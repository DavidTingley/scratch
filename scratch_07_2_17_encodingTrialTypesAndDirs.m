dd  = dir('*201*');

 hpc_corrs =[];
 ls_corrs=[];
  hpc_corrs_shuffle =[];
 ls_corrs_shuffle=[];
 
%% compile data

for ii=15:length(dd)
   cd(dd(ii).name)
   load([dd(ii).name '.behavior.mat'])
   clf
   bz_plotTrials(behavior);
   if ~strcmp(behavior.description,'linear')
           d = input(['directions [l/r]' ': '],'s');
%            tt = input(['types [w/s]' ': '],'s');
       for t = 1:length(unique(behavior.events.trialConditions))
           f = behavior.events.trialConditions(t);
           
           for j=1:length(f)
               if strcmp(d(t),'l')
                   behavior.events.trials{f(j)}.direction = 'counter-clockwise';
               elseif strcmp(d(t),'r')
                   behavior.events.trials{f(j)}.direction = 'clockwise';
               end
           end
               
%            d{t} =  behavior.events.trials{t}.direction;
%            tt{t} =  behavior.events.trials{t}.type;
%            
% %            if strcmp(d(t),'l')
% %                behavior.events.trials{t}.direction = 'counter-clockwise';
% %            elseif strcmp(d(t),'r')
% %                behavior.events.trials{t}.direction = 'clockwise';
% %            end
% %            if strcmp(tt(t),'w')
% %                behavior.events.trials{t}.type = 'wheel alternation';
% %            elseif strcmp(tt(t),'s')
% %                behavior.events.trials{t}.type = 'central alternation';
% %            end
       end
%        for t=1:length(behavior.events.trials)
% %             angles = SpinCalc('QtoEA321',[behavior.events.trials{t}.orientation.x,...
% %             behavior.events.trials{t}.orientation.y,...
% %             behavior.events.trials{t}.orientation.z,...
% %             behavior.events.trials{t}.orientation.w],1,0);
% % %            [p s] = polyfit(length(angles)/2:length(angles),cos(angles(length(angles)/2:end,3))',1)
%            if sum(diff(angles(length(angles)/2:end,3))) < 0
%                behavior.events.trials{t}.direction = 'counter-clockwise';
%            elseif sum(diff(angles(length(angles)/2:end,3))) > 0
%                behavior.events.trials{t}.direction = 'clockwise';
%            end
%            f = behavior.events.trialConditions(t);
% %            behavior.events.trials{t}.direction = d{f};
% %            behavior.events.trials{t}.type = tt{f};
%        end
       save([dd(ii).name '.behavior.mat'],'behavior')
   end
   clear behavior tt d
   cd /home/david/datasets/lsDataset/
end