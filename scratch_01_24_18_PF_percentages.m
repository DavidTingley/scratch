d  = dir('*201*');
for th = 1:5
    hpc_fields{th} = [];
    ls_fields{th} = [];
end
thresh = [1 5 10 20 40];
peakFR = [];
width = [];

for i=1:length(d)
    i
   cd(d(i).name) 
   spikes = bz_GetSpikes('noprompt',true);
   load([d(i).name '.behavior.mat'])
   load([d(i).name '.firingMaps.cellinfo.mat'])
   for th = 4%1:5
      if exist([d(i).name '.placeFields.' num2str(thresh(th),'%0.2d') '_pctThresh.mat'])
      load([d(i).name '.placeFields.' num2str(thresh(th),'%0.2d') '_pctThresh.mat'],'fields')  
      for cond = 1:length(fields)
          if sum(behavior.events.trialConditions==cond) >= 10 %%%%%%%%%%%%%%%%%%%%%%%%%%
          for cell = 1:length(fields{cond})
              if sum(sum(firingMaps.countMaps{cond}(cell,:,:))) >= 1.5 * sum(behavior.events.trialConditions==cond) 
              if ~isempty(fields{cond}{cell})
                  for field = 1:length(fields{cond}{cell})
                      if strcmp(spikes.region{cell},'ls')
                          ls_fields{th} = [ls_fields{th};fields{cond}{cell}{field}.COM];
                      elseif strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca3') | strcmp(spikes.region{cell},'ca1')
                          hpc_fields{th} = [hpc_fields{th};fields{cond}{cell}{field}.COM];
                          peakFR = [peakFR; fields{cond}{cell}{field}.peakFR];
                          width = [width; fields{cond}{cell}{field}.width];
                      end 
                  end
               else
               if strcmp(spikes.region{cell},'ls')
               ls_fields{th} = [ls_fields{th};nan];
               elseif strcmp(spikes.region{cell},'hpc')
               hpc_fields{th} = [hpc_fields{th};nan];
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
   cd /home/david/datasets/lsDataset/
end

for i=1:4
plot(i,mean(~isnan(ls_fields{i})),'.r')
hold on
plot(i,mean(~isnan(hpc_fields{i})),'.k')
end




