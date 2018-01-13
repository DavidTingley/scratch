d = dir('*201*');
count = 1;

for rec = 1:length(d)
   cd(d(rec).name)
   
   sessionInfo = bz_getSessionInfo;
   if exist([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.placeFields.01_pctThresh.mat'])
   spikes = bz_GetSpikes;
   
   for cond = 1:length(phaseMaps.phaseMaps)
      for cell = 1:length(phaseMaps.region) 
          if strcmp(phaseMaps.region{cell},'ls')
          if ~isempty(phaseMaps.phaseMaps{cond}{cell})
            p = phaseMaps.phaseMaps{cond}{cell};
            if size(p,1) > 15
            c=1;
            for i=-pi:.01:pi
                y = bz_wrap(p(:,end),i);
                [a(c) b]=corr(y,p(:,1));
                c=1+c;
            end
            [blah offset] = min(a);
            span = -pi:.01:pi;
            offset = span(offset);
            [a b m] = polyfit(p(:,1),bz_wrap(p(:,end),offset),1);
            slope(count) = a(1);
            region(count) =1;
            additionalDepth = find(sessionInfo.spikeGroups.groups{spikes.shankID(cell)}==spikes.maxWaveformCh(cell))*10;
            depth(count) = [(sessionInfo.depth)+additionalDepth];
            animal(count) = sum(double(sessionInfo.animal));
            recording(count) = rec;
            count = 1 + count;
            end
          end
          elseif strcmp(phaseMaps.region{cell},'hpc') | strcmp(phaseMaps.region{cell},'ca3') |strcmp(phaseMaps.region{cell},'ca1')
          if ~isempty(phaseMaps.phaseMaps{cond}{cell}) & ~isempty(fields{cond}{cell})
            for field = 1:length(fields{cond}{cell})
            p = phaseMaps.phaseMaps{cond}{cell};
            p = Restrict(p,[fields{cond}{cell}{field}.start fields{cond}{cell}{field}.stop]);
            if size(p,1) > 15
            c=1;
            for i=-pi:.01:pi
                y = bz_wrap(p(:,end),i);
                [a(c) b]=corr(y,p(:,1));
                c=1+c;
            end
            [blah offset] = min(a);
            span = -pi:.01:pi;
            offset = span(offset);
            [a b m] = polyfit(p(:,1),bz_wrap(p(:,end),offset),1);
            slope(count) = a(1);
            region(count) = 2;
            additionalDepth = find(sessionInfo.spikeGroups.groups{spikes.shankID(cell)}==spikes.maxWaveformCh(cell))*10;
            depth(count) = [;(sessionInfo.depth)+additionalDepth];
            animal(count) = sum(double(sessionInfo.animal));
            recording(count) = rec;
            count = 1 + count;
            end
            end
          end
          end
      end
   end 
   end
   
   cd('/home/david/datasets/lsDataset')
end