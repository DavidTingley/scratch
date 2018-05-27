d  = dir('*201*');
for th = 1:5
    hpc_fields{th} = [];
    ls_fields{th} = [];
end
thresh = [1 5 10 20 40];

for i=1:length(d)
    i
   cd(d(i).name) 
   spikes = bz_GetSpikes('noprompt',true);
   for th = 1:5
      if exist([d(i).name '.placeFields.' num2str(thresh(th),'%0.2d') '_pctThresh.mat'])
      load([d(i).name '.placeFields.' num2str(thresh(th),'%0.2d') '_pctThresh.mat'],'fields')  
      for cond = 1:length(fields)
          for cell = 1:length(fields{cond})
              if ~isempty(fields{cond}{cell})
                  for field = 1;%1:length(fields{cond}{cell})
                      if strcmp(spikes.region{cell},'ls')
                          ls_fields{th} = [ls_fields{th};fields{cond}{cell}{field}.COM];
                      elseif strcmp(spikes.region{cell},'hpc')
                          hpc_fields{th} = [hpc_fields{th};fields{cond}{cell}{field}.COM];
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




