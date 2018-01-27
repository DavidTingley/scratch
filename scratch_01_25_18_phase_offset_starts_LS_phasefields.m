homeDirectory = pwd;
recordingList = dir('*201*');

%% get rate/phase maps from all sessions
for rec=1:length(recordingList)
   cd(recordingList(rec).name) 
   sessionInfo = bz_getSessionInfo;
   if exist([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.behavior.mat'])
   load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
   spikes = bz_GetSpikes;
   
   conditions = length(unique(behavior.events.trialConditions));
   circ_m = nan(length(spikes.times),conditions,200);
   pval = nan(length(spikes.times),conditions,200);
   
   for cell = 1:length(spikes.times)
   for cond=1:conditions
    for j=1:200
        if ~isempty(phaseMaps.phaseMaps{cond}{cell})
        f = find(phaseMaps.phaseMaps{cond}{cell}(:,1)>j-15);
        ff = find(phaseMaps.phaseMaps{cond}{cell}(:,1)<j+15);
        if ~isempty(intersect(f,ff))
        circ_m(cell,cond,j) = circ_mean(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7));
        [pval(cell,cond,j) z] = circ_rtest(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7));
        end
        end
    end
    % check if is field
    
   jumps = find(diff(pval(cell,cond,:)<.01));
   jumps = unique([1 jumps' 201]);
   phaseField{cell,cond}=[];
   phaseField_stop{cell,cond}=[];
   phaseOffset{cell,cond}=[];
   hasField{cell,cond}=[];
   
   for j=1:length(jumps)-1
       if jumps(j+1) - jumps(j) > 30 & pval(cell,cond,jumps(j)+1) < .01
           if ~isempty(phaseMaps.phaseMaps{cond}{cell})
           phaseField{cell,cond} = [phaseField{cell,cond}; jumps(j)];
           phaseField_stop{cell,cond} = phaseField_stop{cell,cond};jumps(j+1)];
           f = find(phaseMaps.phaseMaps{cond}{cell}(:,1)>jumps(j)-15);
           ff = find(phaseMaps.phaseMaps{cond}{cell}(:,1)<jumps(j)+15);
           phaseOffset{cell,cond} =  [phaseOffset{cell,cond};circ_mean(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7))];
           if ~isempty(fields{cond}{cell})
               hasField{cell,cond} = [hasField{cell,cond};1];
           else
               hasField{cell,cond} = [hasField{cell,cond};0];
           end
           end
       end
   end
   end
   end
   
   offsets{rec} = phaseOffset;
   fieldStarts{rec} = phaseField;
   regions{rec} = phaseMaps.region;
   PF_field{rec} = hasField;
   clear phaseOffset phaseField pval circ_m jumps hasField
   end
   
   cd /home/david/datasets/lsDataset
   
end


%% now lets plot some ish
reg = [];
off = [];
starts = [];
PF = [];

for rec=1:length(recordingList)
    for cell = 1:size(offsets{rec},1)
        for cond=1:size(offsets{rec},2)
            for field = 1:length(offsets{rec}{cell,cond})
                off = [off;offsets{rec}{cell,cond}(field)];
                starts = [starts; fieldStarts{rec}{cell,cond}(field)];
                reg = [reg;sum(double(regions{rec}{cell}))];
                PF = [PF; PF_field{rec}{cell,cond}(field)];
            end
        end
    end
end

subplot(2,2,1)
subplot(2,2,2)
subplot(2,2,3)
subplot(2,2,4)




