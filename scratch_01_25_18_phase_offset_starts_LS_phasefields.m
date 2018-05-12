homeDirectory = pwd;
recordingList = dir('*201*');
winHalfWidth = 10;

%% get rate/phase maps from all sessions
for rec=1:length(recordingList)
   cd(recordingList(rec).name) 
   sessionInfo = bz_getSessionInfo;
   if exist([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.behavior.mat'])
   load([sessionInfo.FileName '.placeFields.20_pctThresh.mat'])
   spikes = bz_GetSpikes;
   load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
   conditions = length(unique(behavior.events.trialConditions));
   circ_m = nan(length(spikes.times),conditions,201);
   pval = nan(length(spikes.times),conditions,201);
   
   for cell = 1:length(spikes.times)
   if strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca3') | strcmp(spikes.region{cell},'ca1')
   for cond=1:conditions
       if ~isempty(fields{cond}{cell})
       row = find(positionDecodingMaxCorr_binned_box_mean.results{cell}.condition == cond);
       col = find(positionDecodingMaxCorr_binned_box_mean.results{cell}.tau == 20);
       row = intersect(row,col);
       [sig_c pval_c] = ttest2(positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_chance_phase(row),positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_phase(row));
    
    for j=1:201
        if ~isempty(phaseMaps.phaseMaps{cond}{cell})
        f = find(phaseMaps.phaseMaps{cond}{cell}(:,1)>j-winHalfWidth);
        ff = find(phaseMaps.phaseMaps{cond}{cell}(:,1)<j+winHalfWidth);
        if ~isempty(intersect(f,ff))
        circ_m(cell,cond,j) = circ_mean(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7));
        [pval(cell,cond,j) z] = circ_rtest(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7));
        end
        end
    end
    % check if is field
    
   jumps = find(diff(pval(cell,cond,:)<.01));
   jumps = unique([1 jumps' 201]);
   phaseDist{cell,cond}=[];
   phaseField{cell,cond}=[];
   phaseField_stop{cell,cond}=[];
   phaseOffset{cell,cond}=[];
   phaseOffset_stop{cell,cond}=[];
   hasField{cell,cond}=[];
   phaseCodingPval{cell,cond} = pval_c;
   
   for j=1:length(jumps)-1
       if jumps(j+1) - jumps(j) > 8 & mean(pval(cell,cond,jumps(j):jumps(j+1)))<.01  & jumps(j+1) - jumps(j) < 100
           if ~isempty(phaseMaps.phaseMaps{cond}{cell})
           phaseField{cell,cond} = [phaseField{cell,cond}; jumps(j)];
           phaseField_stop{cell,cond} = [phaseField_stop{cell,cond};jumps(j+1)];
           f = find(phaseMaps.phaseMaps{cond}{cell}(:,1)>jumps(j)-winHalfWidth);
           ff = find(phaseMaps.phaseMaps{cond}{cell}(:,1)<jumps(j)+winHalfWidth);
           phaseOffset{cell,cond} =  [phaseOffset{cell,cond};circ_mean(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7))];
           
           if all(pval(cell,cond,jumps(j):jumps(j+1))<.05)
               c=1;
               for step =jumps(j):jumps(j+1)
                    f = find(phaseMaps.phaseMaps{cond}{cell}(:,1)>step-winHalfWidth);
                    ff = find(phaseMaps.phaseMaps{cond}{cell}(:,1)<step+winHalfWidth);
                    m(c) = circ_mean(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7));
                    c = 1 + c;          
               end
               subplot(2,2,1)
               plot(unwrap(m))
               hold on
               subplot(2,2,2)
               plot(jumps(j):jumps(j+1),unwrap(m))
               hold on
               subplot(2,2,3)
               plot(unwrap(m)-m(1))
               hold on
%                [p] = polyfit(1:length(m),m,1);
                pause(.01)
               phaseDist{cell,cond} = [phaseDist{cell,cond};sum(diff(unwrap(m)))];
               clear m
           else
               phaseDist{cell,cond} = [phaseDist{cell,cond};nan];
           end
           
           f = find(phaseMaps.phaseMaps{cond}{cell}(:,1)>jumps(j+1)-winHalfWidth);
           ff = find(phaseMaps.phaseMaps{cond}{cell}(:,1)<jumps(j+1)+winHalfWidth);
           phaseOffset_stop{cell,cond} =  [phaseOffset_stop{cell,cond};circ_mean(phaseMaps.phaseMaps{cond}{cell}(intersect(f,ff),7))];
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
   else
       for cond=1:conditions
          phaseDist{cell,cond}=[];
          phaseField{cell,cond}=[];
          phaseField_stop{cell,cond}=[];
          phaseOffset{cell,cond}=[];
          phaseOffset_stop{cell,cond}=[];
          hasField{cell,cond}=[];
          phaseCodingPval{cell,cond} = [];
       end
   end
   end
   
   dist{rec} = phaseDist;
   offsets{rec} = phaseOffset;
   offsets_stop{rec} = phaseOffset_stop;
   fieldStarts{rec} = phaseField;
   fieldStops{rec} = phaseField_stop;
   regions{rec} = phaseMaps.region;
   PF_field{rec} = hasField;
   pv{rec} = phaseCodingPval;
   clear phaseOffset phaseOffset_stop phaseField pval circ_m jumps hasField phaseCodingPval
   end
   
   cd /home/david/datasets/lsDataset
   
end


%% now lets plot some ish
p=[];
reg = [];
distances = [];
off = [];
off_stop = [];
starts = [];
stops = [];
PF = [];

for rec=1:length(recordingList)
    for cell = 1:size(offsets{rec},1)
        for cond=1:size(offsets{rec},2)
            for field = 1:length(offsets{rec}{cell,cond})
                p = [p;pv{rec}{cell,cond}];
                distances =[distances; dist{rec}{cell,cond}(field)];
                off = [off;offsets{rec}{cell,cond}(field)];
                off_stop = [off_stop;offsets_stop{rec}{cell,cond}(field)];
                starts = [starts; fieldStarts{rec}{cell,cond}(field)];
                stops = [stops; fieldStops{rec}{cell,cond}(field)];
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




ls