add

% clf
% 
d  = dir('*201*');
for ii=1:length(d)
   cd(d(ii).name) 
%    load([d(ii).name '.isolationMetrics.cellinfo.mat'])
    if ~exist([d(ii).name '.isolationMetrics.cellinfo.mat'])
        spikes = bz_GetSpikes;
        if ~isempty(spikes)
        isolationMetrics.UID = spikes.UID;
        isolationMetrics.region = spikes.region;
        isolationMetrics.sessionName = spikes.sessionName;

        load([spikes.sessionName '.sessionInfo.mat'])

        for i=1:length(sessionInfo.SpkGrps)
            disp(['working on ' num2str(i)])
            f = find(spikes.shankID==i);
            if exist([sessionInfo.FileName '.fet.' num2str(i)])
                [Fet, nFeatures] = LoadFeatures([sessionInfo.FileName '.fet.' num2str(i)]);
                nChannels = length(sessionInfo.SpkGrps(i).Channels);    
                nSamples = sessionInfo.SpkGrps(i).nSamples;    
                try
                    waveforms = LoadSpikeWaveforms([sessionInfo.FileName '.spk.' num2str(i)],nChannels,nSamples);
                    waveforms = single(reshape(waveforms,length(waveforms),nSamples*nChannels));
                catch
                    waveforms = LoadSpikeWaveforms([sessionInfo.FileName '.spk.' num2str(i)],nChannels*nSamples,1);
                end
                clu = load([[sessionInfo.FileName '.clu.' num2str(i)]]);
                clu = clu(2:end);
                

                for cell=1:length(f)
                    isolationMetrics.isoDist_wav(f(cell)) = IsolationDistance(waveforms,find(clu==spikes.cluID(f(cell))));
                    isolationMetrics.isoDist_fet(f(cell)) = IsolationDistance(Fet,find(clu==spikes.cluID(f(cell))));
                    isolationMetrics.lratio_wav(f(cell)) = L_Ratio(waveforms,find(clu==spikes.cluID(f(cell))));
                    isolationMetrics.lratio_fet(f(cell)) = L_Ratio(Fet,find(clu==spikes.cluID(f(cell))));
                    isolationMetrics.waveformAmplitude(f(cell)) = (min((spikes.rawWaveform{f(cell)}))-max((spikes.rawWaveform{f(cell)}))) *.195;
                end   
            end
        end
        end
        isolationMetrics.waveformUnits = 'uV';
        save([isolationMetrics.sessionName '.isolationMetrics.cellinfo.mat'],'isolationMetrics')
        clear isolationMetrics
    end
 cd /home/david/datasets/lsDataset/
end


% 
% 
% % 
% c=1;
% d = ls('*/*isol*');
% dd = strsplit(d,'\n');
% fdists = [];
% amplitudes = [];
% for i=1:length(dd)-1
%     load(dd{i})
%     isolationMetrics
%     
%     fdists = [fdists,isolationMetrics.isoDist_fet];
%     amplitudes = [amplitudes,isolationMetrics.waveformAmplitude];
%     for j=1:length(isolationMetrics.region)
%         region(c) = sum(double(isolationMetrics.region{j}));
%         c=1+c;
%     end
% end
% 
% subplot(2,2,1)
% histogram(fdists(region==223),0:5:700)
% title(['ls ' num2str(nanmean(fdists(region==223))) ' total:' num2str(length(find(region==223)))])
% set(gca,'xscale','log')
% axis([10 1000 0 350])
% subplot(2,2,2)
% histogram(fdists(region==315),0:5:700)
% title(['hpc ' num2str(nanmean(fdists(region==315))) ' total:' num2str(length(find(region==315)))])
% set(gca,'xscale','log')
% axis([10 1000 0 350])
% 
% 
% subplot(2,2,3)
% histogram(abs(amplitudes(region==223)),0:10:1300)
% title(['ls ' num2str(nanmean(abs(amplitudes(region==223)))) ' total:' num2str(length(find(region==223)))])
% set(gca,'xscale','log')
% axis([10 1000 0 350])
% subplot(2,2,4)
% histogram(abs(amplitudes(region==315)),0:10:1300)
% title(['hpc ' num2str(nanmean(abs(amplitudes(region==315)))) ' total:' num2str(length(find(region==315)))])
% set(gca,'xscale','log')
% axis([10 1000 0 350])
% % 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
