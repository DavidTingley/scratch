clear
clf
cd /home/david/datasets/lsDataset
d  = dir('*201*');
count= 1;

for ii=1:length(d)
cd(d(ii).name) 
    if exist([d(ii).name '.firingMaps.cellinfo.mat'])
        load([d(ii).name '.firingMaps.cellinfo.mat'])
        load([d(ii).name '.placeFields.20_pctThresh.mat'])
        sessionInfo = bz_getSessionInfo;
        conditions = 1:length(firingMaps.rateMaps);
        for cell =1:size(firingMaps.rateMaps{1},1)
        for cond = 1:conditions
            meanRate(count,:) = squeeze(mean(firingMaps.rateMaps{cond}(cell,:,:)));
%         for field = 1:length(fields{cond}{cell})
%             ph = firingMaps.phaseMaps{cond}{cell}(:,[1 7]);
%             start = fields{cond}{cell}{field}.start;
%             stop = fields{cond}{cell}{field}.stop;
%             p = Restrict(ph,[start stop]);
%             if size(p,1) > 15
%                 [P,S] = polyfit(p(:,1),cos(p(:,2)),1);
%                 slope(count) = P(1);
%                 correlation(count) = circ_corrcl(p(:,2),p(:,1));
%                 maxR(count) = max(firingMaps.phaseMaps{cond}{cell}(:,[6]));
%                 com(count) = fields{cond}{cell}{field}.COM;
                if isempty(sessionInfo.ca3)
                    region(count) = 1;
                else
                    region(count) = 3;
                end
%                 subplot(2,2,region(count))
%                 scatter(maxR(count),correlation(count),'.k')
%                 hold on
%                 subplot(2,2,region(count)+1)
%                 scatter(com(count),maxR(count),'.r')
%                 hold on
%                 pause(.01)
                count = 1+count;
%             end
%         end
        end
        end
    end
cd /home/david/datasets/lsDataset
end
