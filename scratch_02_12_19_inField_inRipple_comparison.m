clear all
cd /home/david/datasets/ripples_LS/
d = dir('*201*');
nCells = 1;

% can the same analysis be done with population FR? could be a more
% accurate read out at the population level..

for i=1:length(d)
   cd(d(i).name)
   
%    ripples = bz_LoadEvents(pwd,'CA1Ripples');
   ripples = bz_LoadEvents(pwd,'popBursts');
   if exist([d(i).name '.firingMaps.cellinfo.mat'])
   maps = load([d(i).name '.firingMaps.cellinfo.mat']);
   fields = load([d(i).name '.placeFields.20_pctThresh.mat']);
   spikes = bz_GetSpikes('noprompts',true);
   
   if ~isempty(ripples) & ~isempty(maps) & ~isempty(fields)
       
       
       for sp = 1:length(spikes.times)
           if strcmp(spikes.region{sp},'hpc') | strcmp(spikes.region{sp},'ca3') | strcmp(spikes.region{sp},'ca1')
               spkFieldCount{nCells} = [];
               spkFieldRate{nCells} = [];
               spkRipCount{nCells} = [];
               spkRipRate{nCells} = [];
       
               for cond = 1:length(fields.fields)
                   if ~isempty(fields.fields{cond}{sp})
                       start = fields.fields{cond}{sp}{1}.start;
                       stop =  fields.fields{cond}{sp}{1}.stop;
                       for trial = 1:size(maps.firingMaps.countMaps{cond},2)
%                           spkFieldCount{sp}{cond}(trial) = sum(maps.firingMaps.countMaps{cond}(sp,trial,start:stop));
                          spkFieldRate{nCells}{cond}(trial) = mean(maps.firingMaps.rateMaps{cond}(sp,trial,start:stop));
                          spkFieldCount{nCells}{cond}(trial) = sum(maps.firingMaps.countMaps{cond}(sp,trial,start:stop));
                       end
                   end
               end
               for rip =1:length(ripples.timestamps)
                   spks = Restrict(spikes.times{sp},[ripples.timestamps(rip,1) ripples.timestamps(rip,2)]);
%                    spkRipCount{sp}(rip) = length(spks);
                   spkRipRate{nCells}(rip) = length(spks)./(diff(ripples.timestamps(rip,:)));
                   spkRipCount{nCells}(rip) = length(spks);
               end
               
            subplot(2,2,1)
            scatter(mean(spkRipCount{nCells}(spkRipCount{nCells}>0)),mean(cell2mat(spkFieldCount{nCells})),'.k')
            hold on
            subplot(2,2,2)
            scatter(mean(spkRipRate{nCells}(spkRipRate{nCells}>0)),mean(cell2mat(spkFieldRate{nCells})),'.k')
            hold on
            subplot(2,2,3)

            subplot(2,2,4)  
            pause(.01)
            
            nCells = 1+nCells;
            
           end
       end

   
   end
   end
   

   
   
   cd /home/david/datasets/ripples_LS
    
    
    
end