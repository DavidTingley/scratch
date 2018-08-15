
cd D:\Dropbox\datasets\lsDataset
ls_fields = [];
hpc_fields = []; 
peakFR = [];
width = [];
c=1;
d  = dir('*201*');
thresh = 05;

for i=1:length(d)
    
   cd(d(i).name) 
   if exist([d(i).name '.firingMaps.cellinfo.mat'])
   spikes = bz_GetSpikes('noprompt',true);
   load([d(i).name '.behavior.mat'])
   load([d(i).name '.firingMaps.cellinfo.mat'])
      if exist([d(i).name '.placeFields.' num2str(thresh,'%0.2d') '_pctThresh.mat'])
      load([d(i).name '.placeFields.' num2str(thresh,'%0.2d') '_pctThresh.mat'],'fields')  
      % get speed
      for trial = 1:length(behavior.events.trials)
          behavior.events.trials{trial}.speed = [smooth(abs(diff(behavior.events.trials{trial}.x))+abs(diff(behavior.events.trials{trial}.y)),15);nan];
      end
       
       
      for cond = 1:length(fields)
          if sum(behavior.events.trialConditions==cond) >= 10 && ...
                  strcmp(behavior.events.conditionType{cond},'wheel') %%%%%%%%%%%%%%%%%%%%%%%%%%
          for cell = 1:length(fields{cond})
              if sum(sum(firingMaps.countMaps{cond}(cell,:,:))) >= 1.5 * sum(behavior.events.trialConditions==cond) 
              if ~isempty(fields{cond}{cell})
                  for field = 1:length(fields{cond}{cell})
                      if strcmp(spikes.region{cell},'ls')
                          ls_fields = [ls_fields;fields{cond}{cell}{field}.COM];
                      elseif strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca3') | strcmp(spikes.region{cell},'ca1')
                          trials = find(behavior.events.trialConditions==cond);
                              hpc_fields = [hpc_fields;fields{cond}{cell}{field}.COM];
                              peakFR = [peakFR; fields{cond}{cell}{field}.peakFR];
                              width = [width; fields{cond}{cell}{field}.width];


                              for trial = 1:length(trials)
                                  start = find(behavior.events.trials{trials(trial)}.mapping == fields{cond}{cell}{field}.start,1,'first');
                                  stop = find(behavior.events.trials{trials(trial)}.mapping == fields{cond}{cell}{field}.stop,1,'last');
                                  if isempty(start) | isempty(stop)
                                    [nah start] = min(abs(behavior.events.trials{trials(trial)}.mapping - fields{cond}{cell}{field}.start));
                                    [nah stop] = min(abs(behavior.events.trials{trials(trial)}.mapping - fields{cond}{cell}{field}.stop));
                                  end
                                  speeds(trial) = nanmean(behavior.events.trials{trials(trial)}.speed(start:stop));

                                  [nah peaks(trial)] = max(firingMaps.rateMaps{cond}(cell,trial,fields{cond}{cell}{field}.start:fields{cond}{cell}{field}.stop));
                                  peaks(trial) = peaks(trial) + fields{cond}{cell}{field}.start;

                              end
                              if sum(~isnan(speeds)) ~= length(trials)
                                  error
                              end
                              scatter(peaks,speeds,'.k')
                              axis([0 200 0 max(speeds)+.000001])
                              peaks(peaks==0)=nan; peaks(peaks==length(peaks))=nan;
                              [co(c) pval(c)] = corr(peaks',speeds');
                              c=1+c;
                              pause(.1); clear peaks speeds
                              
                      end 
                  end
               else
               if strcmp(spikes.region{cell},'ls')
               ls_fields = [ls_fields;nan];
               elseif strcmp(spikes.region{cell},'hpc')
               hpc_fields = [hpc_fields;nan];
               co(c) =nan;
               pval(c)=nan;
               c=1+c;
               end
              end 
              end
           end
           end
      end
%       else
%           error('couldnt find PF file')
      end
   end
   cd D:\Dropbox\datasets\lsDataset
end






