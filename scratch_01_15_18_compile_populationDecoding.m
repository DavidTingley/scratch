d = dir('*201*');
% clf
tau = 25;
for nCells = 10
%     figure(nCells)
LS=[];HPC=[];REC=[];

for i=1:length(d)
    cd(d(i).name)
    if length(dir('*popinfo*'))==3
    sessionInfo = bz_getSessionInfo;
    load([sessionInfo.FileName '.behavior.mat']);
    ls = load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_LS.popinfo.mat']);
    hpc = load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_HPC.popinfo.mat']);
    rec = load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box.popinfo.mat']);

    conditions = unique(behavior.events.trialConditions);
    for cond = conditions
        if length(behavior.events.trialConditions==cond)>=10
        rows = find(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.tau == tau);
        cols = find(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.condition == cond);
        ind = intersect(rows,cols);
        ls_mse_phase(1) = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase(ind));
        ls_mse_phase(2) = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_all(ind));
        ls_mse_phase(3) = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_cos(ind));
        ls_mse_phase = min(ls_mse_phase);
        ls_mse_rate = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_rate(ind));
        ls_nCells = length(ls.positionDecodingMaxCorr_binned_box_mean.region);
        
        rows = find(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.tau == tau);
        cols = find(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.condition == cond);
        ind = intersect(rows,cols);
        hpc_mse_phase(1) = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase(ind));
        hpc_mse_phase(2) = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_all(ind));
        hpc_mse_phase(3) = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_cos(ind));
        hpc_mse_phase = min(hpc_mse_phase); 
        hpc_mse_rate = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_rate(ind));
        hpc_nCells = length(hpc.positionDecodingMaxCorr_binned_box_mean.region);
        
        rows = find(rec.positionDecodingMaxCorr_binned_box_mean.results{end}.tau == tau);
        cols = find(rec.positionDecodingMaxCorr_binned_box_mean.results{end}.condition == cond);
        ind = intersect(rows,cols);
        rec_mse_phase(1) = nanmedian(rec.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase(ind));
        rec_mse_phase(2) = nanmedian(rec.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_all(ind));
        rec_mse_phase(3) = nanmedian(rec.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_cos(ind));
        rec_mse_phase = min(rec_mse_phase);         
        rec_mse_rate = nanmedian(rec.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_rate(ind));
        rec_nCells = length(rec.positionDecodingMaxCorr_binned_box_mean.region);
    
        subplot(3,3,1)
        scatter(ls_nCells,ls_mse_phase,'.b')
        hold on
        scatter(ls_nCells,ls_mse_rate,'.r')
        axis([0 200 -1000 8000])
        title('ls cells')
        LS = [LS; ls_nCells, ls_mse_phase, ls_mse_rate, i, cond];
        
        subplot(3,3,2)
        scatter(hpc_nCells,hpc_mse_phase,'.b')
        hold on
        scatter(hpc_nCells,hpc_mse_rate,'.r')
        axis([0 200 -1000 8000])
        title('hpc cells')
        HPC = [HPC; hpc_nCells, hpc_mse_phase, hpc_mse_rate, i, cond];
        
        subplot(3,3,3)
        scatter(rec_nCells,rec_mse_phase,'.b')
        hold on
        scatter(rec_nCells,rec_mse_rate,'.r')
        axis([0 200 -1000 8000])
        title('all cells')
        REC = [REC; rec_nCells, rec_mse_phase, rec_mse_rate,  i, cond];
        
        pause(.01)
        end
    end
    end
    cd('D:\Dropbox\datasets\lsDataset')
end

LS(LS(:,1)<nCells,:)=[];



LS(isnan(LS(:,2)),:) = [];
HPC(isnan(HPC(:,2)),:) = [];
REC(isnan(REC(:,2)),:) = [];

for p = 1:2
    subplot(3,3,1)
    [x y] = polyfit(LS(:,1),LS(:,2),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'b')
    hold on
    [x y] = polyfit(LS(:,1),LS(:,3),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'r')
    axis([0 200 -1000 8000])

    subplot(3,3,2)
    [x y] = polyfit(HPC(:,1),HPC(:,2),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'b')
    hold on
    [x y] = polyfit(HPC(:,1),HPC(:,3),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'r')
    axis([0 200 -1000 8000])

    subplot(3,3,3)
    [x y] = polyfit(REC(:,1),REC(:,2),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'b')
    hold on
    [x y] = polyfit(REC(:,1),REC(:,3),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'r')
    axis([0 200 -1000 8000])
end

for i=1:length(d)
    f = find(LS(:,4)==i);
    l(i,:)=nanmedian(LS(f,1:3));
    f = find(LS(:,4)==i);
    h(i,:)=nanmedian(HPC(f,1:3));
    f = find(LS(:,4)==i);
    r(i,:)=nanmedian(REC(f,1:3));
end
LS = l;
HPC=h;
REC=r;

subplot(3,3,7)
scatter(LS(:,1),LS(:,2),'.b')
hold on
scatter(LS(:,1),LS(:,3),'.r')
axis([0 200 -1000 8000])

subplot(3,3,8)
scatter(HPC(:,1),HPC(:,2),'.b')
hold on
scatter(HPC(:,1),LS(:,3),'.r')
axis([0 200 -1000 8000])

subplot(3,3,9)
scatter(REC(:,1),REC(:,2),'.b')
hold on
scatter(REC(:,1),REC(:,3),'.r')
axis([0 200 -1000 8000])
% [x y] = polyfit(depths(f),(lp(f))-(lr(f)),1);
% y1 = polyval(x,depths(f));

% LS(isnan(LS(:,2)),:) = [];
% HPC(isnan(HPC(:,2)),:) = [];
% REC(isnan(REC(:,2)),:) = [];

for p = 1:2
    subplot(3,3,7)
    [x y] = polyfit(LS(:,1),LS(:,2),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'b')
    hold on
    [x y] = polyfit(LS(:,1),LS(:,3),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'r')
    axis([0 200 -1000 8000])

    subplot(3,3,8)
    [x y] = polyfit(HPC(:,1),HPC(:,2),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'b')
    hold on
    [x y] = polyfit(HPC(:,1),HPC(:,3),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'r')
    axis([0 200 -1000 8000])

    subplot(3,3,9)
    [x y] = polyfit(REC(:,1),REC(:,2),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'b')
    hold on
    [x y] = polyfit(REC(:,1),REC(:,3),p);
    y1 = polyval(x,1:200);
    plot(1:200,y1,'r')
    axis([0 200 -1000 8000])
end

subplot(3,3,5)
[x y] = polyfit(LS(:,1),LS(:,2),1);
y1 = polyval(x,1:200);
plot(1:200,y1,'.b')
hold on
[x y] = polyfit(LS(:,1),LS(:,3),1);
y1 = polyval(x,1:200);
plot(1:200,y1,'.r')
axis([0 200 -1000 8000])
[x y] = polyfit(HPC(:,1),HPC(:,2),1);
y1 = polyval(x,1:200);
plot(1:200,y1,'b')
hold on
[x y] = polyfit(HPC(:,1),HPC(:,3),1);
y1 = polyval(x,1:200);
plot(1:200,y1,'r')
axis([0 200 -1000 8000])
[x y] = polyfit(REC(:,1),REC(:,2),1);
y1 = polyval(x,1:200);
plot(1:200,y1,'.k')
hold on
[x y] = polyfit(REC(:,1),REC(:,3),1);
y1 = polyval(x,1:200);
plot(1:200,y1,'k')
axis([0 200 -1000 8000])
end