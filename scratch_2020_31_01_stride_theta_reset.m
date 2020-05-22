folders = dir('*201*');

c=1;
        
parfor f =1:55
    [mats(f)] = runPar(folders,f);
end


function [mats] = runPar(folders,f)

    f
    cd(folders(f).name)
    cc = 1;
    sessionInfo = bz_getSessionInfo;
    mats.test = [];
    if ~isempty(sessionInfo.ca1)
        load([sessionInfo.FileName '.behavior.mat'])

        lfp = bz_GetLFP(sessionInfo.ca1,'intervals',behavior.events.trialIntervals);
        [b a] = butter(3,[5/625 12/625],'bandpass');
%         filt = FiltFiltM(b,a,single(lfp.data));

        for i=1:length(behavior.events.trialIntervals)
            
%             int = behavior.events.trialIntervals(i,:);
            m1 = makelength(diff(behavior.events.trials{i}.z),800);
%             m1 = Smooth(m1,3);
            m2 = makelength(single(lfp(i).data),800);
            m2_f = makelength(FiltFiltM(b,a,double(lfp(i).data)),800);
            mats.mat1(cc,:) = m1;
            mats.mat2(cc,:) = m2;
            mats.mat2_f(cc,:) = m2_f;
            [pks locs] =findpeaks(diff(behavior.events.trials{i}.z),'MinPeakHeight',.45,'MinPeakDistance',10);
            for l = 1:length(locs)
                int1 =  behavior.events.trials{i}.timestamps(locs(l));
                pkTrigThet(c,:) = lfp.data(round(int1*1250-625):round(int1*1250+625));
                pkTrigPhase(c,:) = angle(hilbert(pkTrigThet(c,:)));
                tType{c} = behavior.events.conditionType{behavior.events.trialConditions(i)};
                c=1+c;
            end
            cc=cc+1;
        end
    end
    
    cd /home/david/datasets/lsDataset/ %D:\datasets\lsDataset %
return 
end



% m1=[];
% m=[];
% for i=1:40
%     m=[m;mats(i).mat2_f];
%     m1 = [m1; mats(i).mat1];
% end

% clear c
% parfor i=1:size(mat1,1)
% c(i,:) = ccgBinned(mat1(i,:),mat2(i,:),200);
% end

% save('/home/david/Dropbox/data/maze_stride_theta_reset.mat','-v7.3')
% 
% ints = behavior.events.trialIntervals;
% 
% for cell =17:148
%     subplot(2,2,1)
%     imagesc(squeeze(firingMaps.rateMaps{1}(cell,:,:)))
%     c=1;
%     subplot(2,2,3)
%     cla
%     ph = [];
%     for i=1:160
%         zd = [(angle(hilbert(Smooth(zscore(diff(behavior.events.trials{i}.z)),2))));0];
%         bts = (behavior.events.trials{i}.timestamps);
%         spks = find(InIntervals(spikes.times{cell},ints(i,:)));
%         for s = 1:length(spks)
%             [a b] = min(abs(bts-spikes.times{cell}(spks(s))));
%             scatter(b,zd(b)+2*pi,'.k');
%             hold on
%             scatter(b,zd(b),'.k');
%             h(c) = zd(b);
%             c=c+1;
%         end
%         ph = [ph;zd];
%     end
%     xlim([0 200])
%     subplot(2,2,4)
%     if ~isempty(phaseMaps.phaseMaps{1}{cell})
%         scatter(phaseMaps.phaseMaps{1}{cell}(:,1),phaseMaps.phaseMaps{1}{cell}(:,end),'.k')
%         hold on
%         scatter(phaseMaps.phaseMaps{1}{cell}(:,1),phaseMaps.phaseMaps{1}{cell}(:,end)+2*pi,'.k')
%         xlim([0 200])
%     end
%     hold off
%     if exist('h')
%         occH = hist(ph(:),-pi:.5:pi);
%         spkH = hist(h,-pi:.5:pi); clear h ph
%         subplot(2,2,2)
%         plot(spkH./occH)
%         pause
%     end
% end