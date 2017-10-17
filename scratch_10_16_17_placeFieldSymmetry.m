clear all
d = dir('*/*placeField*');
f = dir('*/*firingMaps*');
c=1;

for i =1:length(d)
    cd(d(i).folder)
   load([d(i).folder '/' d(i).name]) 
   sessionInfo = bz_getSessionInfo;
   load([d(i).folder '/' sessionInfo.FileName '.firingMaps.cellinfo.mat'])
   for cell = 1:length(firingMaps.region)
       for cond = 1:length(firingMaps.rateMaps)
%            fields{cond} = bz_getPlaceFields1D(firingMaps.rateMaps{cond},'minpeakrate',2,'percentThreshold',.1);
           
           if size(fields{cond},2)>=cell
           if ~strcmp(firingMaps.region{cell},'ls') & ~isempty(fields{cond}{cell})
           for field = 1:length(fields{cond}{cell})
               if fields{cond}{cell}{field}.start > 5 & fields{cond}{cell}{field}.stop < 195
               if size(firingMaps.rateMaps{cond},2) > 10
               meanRate = squeeze(mean(firingMaps.rateMaps{cond}(cell,:,:)));
                  plot(meanRate)
                  if mean(meanRate)>100
                     disp() 
                  end
                  com = fields{cond}{cell}{field}.COM;
                  start = fields{cond}{cell}{field}.start;
                  stop = fields{cond}{cell}{field}.stop;
                  frai(c) = (mean(meanRate(start:com)) - mean(meanRate(com+1:stop))) / (mean(meanRate(start:com)) + mean(meanRate(com+1:stop)));
                  skew(c) = skewness(meanRate(start:stop));
                  sk(c) = mean(meanRate(start:stop)-mean(meanRate(start:stop)).^3) / std(meanRate(start:stop)).^3;
                  r(c,:) = meanRate;
                  COM(c) = com;
                  title([num2str(frai(c)) '     '   num2str(sk(c))])
                  pause(.1)
                  c=1+c;
               end
               end
           end
           end
           end
       end
   end
    clear fields firingMaps
    cd /home/david/datasets/lsDataset
end