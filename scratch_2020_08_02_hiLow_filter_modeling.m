warning off
clear all
cgm_files = dir('_*mat');
hourHz = 1/(60*60);
dayHz = 1/(24*60*60);
nyquist = 5./(24*60)/2;

% % figure(1)
for ff=1:length(cgm_files)
    data{ff} = load(cgm_files(ff).name,'idx','zt','theta*','amps','freq','dur','count','isa','isig_levels','spSlope','states','emgSig','mov');
%     temp = load(cgm_files(ff).name(2:end),'idx');
%     data{ff}.idx = temp.idx;
end


%% 80/20 split testing dataa
for iter = 1:1000
    for ff = 1:length(data) % animals
        if isfield(data{ff},'count')
        count = data{ff}.count; 
        isig_levels = data{ff}.isig_levels;
        theta_resid = data{ff}.theta_resid; 
        theta_z = data{ff}.theta_z;
        theta_deltaR = data{ff}.theta_deltaR;
        spSlope = data{ff}.spSlope;
        emgSig = data{ff}.emgSig;
        states = data{ff}.states;
        zeitTime = discretize(data{ff}.zt,24);
        amps = data{ff}.amps;
        dur = data{ff}.dur;
        freq = data{ff}.freq;
        
        if isfield(data{ff},'mov')
            movement = data{ff}.mov;
            bad_movement = 0;
        else
            movement =  rand(1,length(count));%ones(1,length(count));
            bad_movement = 1;
        end
        
%         pred = [nanZscore(theta_resid'),nanZscore(theta_z'),nanZscore(theta_deltaR'),nanZscore(count'),nanZscore(spSlope'),nanZscore(states'),nanZscore(movement'),nanZscore(emgSig'),nanZscore(zeitTime')]; %,nanZscore(zeitTime')   ,nanZscore(theta_z'),amps',dur',freq',;
%         names = [{'theta_resid'},{'theta_z'},{'theta_deltaR'},{'count'},{'spSlope'},{'states'},{'movement'},{'emgSig'},{'ZeitTime'}]; % ,{'ZeitTime'},    {'theta_z'},'amps','dur','freq',
        pred = [nanZscore(theta_resid'),nanZscore(theta_deltaR'),nanZscore(count'),nanZscore(spSlope'),nanZscore(states'),nanZscore(movement'),nanZscore(emgSig')]; %,nanZscore(zeitTime')   ,nanZscore(theta_z'),amps',dur',freq',;
        names = [{'theta_resid'},{'theta_deltaR'},{'count'},{'spSlope'},{'states'},{'movement'},{'emgSig'}]; % ,{'ZeitTime'},    {'theta_z'},'amps','dur','freq',


%         ica = fast_ica(pred(good,:),5,100);
%         [coeff score lat ts explained] = pca(pred);
%         
%         good = find(~isnan(sum(pred')'));
%         bad = find(isnan(sum(pred')'));
%         ica = fast_ica(score(good,1:5),5);
% 
%         pred(good,end+1:end+3) = ica(:,1:3);
        
%         pred(bad,:) = nan;
%         pred(:,end+1:end+2) = score(:,1:2);
        
        for i=1:size(pred(:,1))
            if  ~ismember(i,data{ff}.idx) | isnan(sum(pred(i,:)))
                pred(i,:) = nan;
            end
%             if  ~ismember(i,data{ff}.idx) 
%                 pred(i,:) = nan;
%             end
        end

    trainSplit = randperm(length(count));
    testSplit = trainSplit(1:round(length(count)/3));
    trainSplit(1:round(length(count)/3)) = [];
%         isa = fillmissing((data{ff}.isig_levels),'linear');
%         % isa2 = highpass(isa,dayHz,nyquist*2);
%         hourHz = 1/(60*60);
%         nyquist = 5./(24*60)/2;
%         [b a] = butter(4,[ hourHz./3 / nyquist],'high');
%         isa = filtfilt(b,a,isa);
%         isa(isnan(isig_levels))=nan;
    isa = [0 diff(data{ff}.isig_levels)];
    % glucFlux = highpass(isa,hourHz,nyadquist*2);
    glucFlux = nanZscore(isa);
    for i=2 % temporal offset (5 min bins)
   
        gluc = circshift(glucFlux,-i);
        glucTrain = gluc(trainSplit);
        glucTest = gluc(testSplit);
        
        for ii=1:500
            p = count-movmean(count,ii,'omitnan');
            mdl = fitlm(p(trainSplit)',glucTrain);
            yfit = predict(mdl,p(testSplit)');
            mse_hipass(ff,ii,iter) = nanmean((glucTest-yfit').^2);
            
            pp = medfilt1(count,ii,'omitnan');
            mdl = fitlm(pp(trainSplit)',glucTrain);
            yfit = predict(mdl,pp(testSplit)');
            mse_lowpass(ff,ii,iter) = nanmean((glucTest-yfit').^2);
        end
        subplot(2,2,1)
        imagesc(zscore(squeeze(nanmean(mse_hipass,3)),[],2))
        subplot(2,2,2)
        imagesc(zscore(squeeze(nanmean(mse_lowpass,3)),[],2))
        subplot(2,2,3)
        plot(mean(zscore(squeeze(nanmean(mse_hipass,3)),[],2)))
        subplot(2,2,4)
        plot(mean(zscore(squeeze(nanmean(mse_lowpass,3)),[],2)))
        pause(.01);
    end
        end
    end
end
