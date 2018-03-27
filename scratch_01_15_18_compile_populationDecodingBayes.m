clear 
d = dir('*201*');
clf
LS=[];HPC=[];
% converge faster with random sampling
% ord = randperm(length(d));
LS_mat_phase = nan(150,200);
LS_mat_rate = nan(150,200);
HPC_mat_phase = nan(150,200);
HPC_mat_rate = nan(150,200);
wind = 20; % smoothing window to plot
nCellThresh = 0; % min num cells to be counted
ord = randperm(length(d));

for i=length(d):-1:1
%     cd(d(i).name)
    load(d(ord(i)).name)
for tau = 20
    hpc=[];ls=[];
    if strcmp(posDecodeBayes.region{1},'ls')
        ls = load(d(ord(i)).name);
        conditions = unique(ls.posDecodeBayes.results{end}.condition);
    elseif strcmp(posDecodeBayes.region{1},'hpc')
        hpc = load(d(ord(i)).name);
        conditions = unique(hpc.posDecodeBayes.results{end}.condition);
    end

    for cond = conditions'
%         if length(behavior.events.trialConditions==cond)>=10
        if ~isempty(ls)
        rows = find(ls.posDecodeBayes.results{end}.tau == tau);
        cols = find(ls.posDecodeBayes.results{end}.condition == cond);
        ind = intersect(rows,cols);
        ls_mse_phase(1) = nanmedian(ls.posDecodeBayes.results{end}.mse_phase(ind));
        ls_mse_phase(2) = nanmedian(ls.posDecodeBayes.results{end}.mse_phase_all(ind));
        ls_mse_phase(3) = nanmedian(ls.posDecodeBayes.results{end}.mse_phase_cos(ind));
        ls_mse_phase = min(ls_mse_phase);
        ls_mse_rate = nanmedian(ls.posDecodeBayes.results{end}.mse_rate(ind));
        ls_nCells = length(ls.posDecodeBayes.region);
        elseif ~isempty(hpc)
        rows = find(hpc.posDecodeBayes.results{end}.tau == tau);
        cols = find(hpc.posDecodeBayes.results{end}.condition == cond);
        ind = intersect(rows,cols);
        hpc_mse_phase(1) = nanmedian(hpc.posDecodeBayes.results{end}.mse_phase(ind));
        hpc_mse_phase(2) = nanmedian(hpc.posDecodeBayes.results{end}.mse_phase_all(ind));
        hpc_mse_phase(3) = nanmedian(hpc.posDecodeBayes.results{end}.mse_phase_cos(ind));
        hpc_mse_phase = min(hpc_mse_phase); 
        hpc_mse_rate = nanmedian(hpc.posDecodeBayes.results{end}.mse_rate(ind));
        hpc_nCells = length(hpc.posDecodeBayes.region);
        end

%         if ~isempty(ls) & tau == wind
%                 subplot(3,3,1)
%                 scatter(ls_nCells,ls_mse_phase,'.b')
%                 hold on
%                 scatter(ls_nCells,ls_mse_rate,'.r')
%                 axis([0 100 -100 14000])
%                 title('ls cells')
%         end
%         if ~isempty(hpc) & tau == wind  
%                 subplot(3,3,2)
%                 scatter(hpc_nCells,hpc_mse_phase,'.b')
%                 hold on
%                 scatter(hpc_nCells,hpc_mse_rate,'.r')
%                 axis([0 100 -100 14000])
%                 title('hpc cells')
%         end
        if ~isempty(ls) 
                LS = [LS; ls_nCells, tau, ls_mse_phase, ls_mse_rate];
        end
        if ~isempty(hpc)
                HPC = [HPC; hpc_nCells, tau, hpc_mse_phase, hpc_mse_rate];
        end
        pause(.01)
        
    end
end
    i
    if ~mod(i,5000)
        if ~isempty(LS)
        LS(isnan(LS(:,3)),:) = [];
        end
        if ~isempty(HPC)
        HPC(isnan(HPC(:,3)),:) = [];
        end
        % REC(isnan(REC(:,2)),:) = [];
%         subplot(3,3,4)
%         hold off
%         subplot(3,3,5)
%         hold off
%         method = {'poly1','exp1','exp2'};
%         for p = 1:2
%             if ~isempty(LS)
%             ind = LS(:,2)==wind & LS(:,1) > nCellThresh;
%             if sum(ind) > 2
%             subplot(3,3,4)
%             [f] = fit(LS(ind,1),LS(ind,3),method{p});
% %             y1 = polyval(x,1:200);
%             plot(f(1:200),'b')
%             hold on
%             [f] = fit(LS(ind,1),LS(ind,4),method{p});
% %             y1 = polyval(x,1:200);
%             plot(f(1:200),'r')
%             axis([0 150 -100 8000])
%             end
%             end
%             if ~isempty(HPC)
%             ind = HPC(:,2)==wind & HPC(:,1) > nCellThresh;
%             if sum(ind)>2
%             subplot(3,3,5)
%             [f] = fit(HPC(ind,1),HPC(ind,3),method{p});
% %             y1 = polyval(x,1:200);
%             plot(f(1:200),'b')
%             hold on
%             [f] = fit(HPC(ind,1),HPC(ind,4),method{p});
% %             y1 = polyval(x,1:200);
%             plot(f(1:200),'r')
%             axis([0 150 -100 8000])
%             end
%             end
%         end
        if ~isempty(LS)
        LS_mat_phase = griddata(LS(:,1),LS(:,2),LS(:,3),meshgrid(1:200,1:200),meshgrid(1:200,1:200));
%         LS_mat_phase = f(meshgrid(1:200,1:200),meshgrid(1:200,1:200));
        LS_mat_phase(LS_mat_phase<0)=0;
        LS_mat_rate = griddata(LS(:,1),LS(:,2),LS(:,4),meshgrid(1:200,1:200),meshgrid(1:200,1:200));
%         LS_mat_rate = f(meshgrid(1:200,1:200),meshgrid(1:200,1:200));
        LS_mat_rate(LS_mat_rate<0)=0;
        end
        if ~isempty(HPC)
        HPC_mat_phase = griddata(HPC(:,1),HPC(:,2),HPC(:,3),meshgrid(1:200,1:200),meshgrid(1:200,1:200));
%         HPC_mat_phase = f(meshgrid(1:200,1:200),meshgrid(1:200,1:200));
        HPC_mat_phase(HPC_mat_phase<0)=0;
        HPC_mat_rate = griddata(HPC(:,1),HPC(:,2),HPC(:,4),meshgrid(1:200,1:200),meshgrid(1:200,1:200));
%         HPC_mat_rate = f(meshgrid(1:200,1:200),meshgrid(1:200,1:200));
        HPC_mat_rate(HPC_mat_rate<0)=0;
        end

        for ii=1:150 % cell num
            if ~isempty(LS)
            f = find(LS(:,1)==ii);
            end
            if ~isempty(HPC)
                g = find(HPC(:,1)==ii);
            end
            for j=1:200 % tau
                % set minimum number of models to be included?
                if ~isempty(LS)
                ff = find(LS(:,2)==j);
                LS_mat_phase(ii,j) = nanmedian(LS(intersect(f,ff),3));
                LS_mat_rate(ii,j) = nanmedian(LS(intersect(f,ff),4));
                LS_mat_phase_ste(ii,j) = nanstd(LS(intersect(f,ff),3))./sqrt(length(intersect(f,ff)));
                LS_mat_rate_ste(ii,j) = nanstd(LS(intersect(f,ff),4))./sqrt(length(intersect(f,ff)));
                end
                if ~isempty(HPC)
                gg = find(HPC(:,2)==j);
                HPC_mat_phase(ii,j) = nanmedian(HPC(intersect(g,gg),3));
                HPC_mat_rate(ii,j) = nanmedian(HPC(intersect(g,gg),4));
                HPC_mat_phase_ste(ii,j) = nanstd(HPC(intersect(g,gg),3))./sqrt(length(intersect(g,gg)));
                HPC_mat_rate_ste(ii,j) = nanstd(HPC(intersect(g,gg),4))./sqrt(length(intersect(g,gg)));
                end
            end
%             LS_mat_phase(ii,:)=fillmissing(LS_mat_phase(ii,:),'linear');
%             LS_mat_rate(ii,:)=fillmissing(LS_mat_rate(ii,:),'linear');
%             if ~isempty(HPC)
%             HPC_mat_phase(ii,:)=fillmissing(HPC_mat_phase(ii,:),'linear');
%             HPC_mat_rate(ii,:)=fillmissing(HPC_mat_rate(ii,:),'linear');
%             end
        end
        if ~isempty(LS)
%         LS_mat_phase(LS_mat_phase<0)=0;
%         LS_mat_rate(LS_mat_rate<0)=0;
        subplot(3,3,6)
        imagesc(log(LS_mat_phase))
%         caxis([0 log(8000)])
        ylabel('cell #')
        xlabel('smoothing window')
        subplot(3,3,7)
        imagesc(log(LS_mat_rate))
%         caxis([0 log(8000)])
        subplot(3,3,8)
        hold off
        plot(1)
%         boundedline(1:150,fillmissing(LS_mat_phase(:,wind),'nearest'),fillmissing(LS_mat_phase_ste(:,wind),'nearest'),'b','transparency',.1)
%         boundedline(1:150,fillmissing(LS_mat_rate(:,wind),'nearest'),fillmissing(LS_mat_rate_ste(:,wind),'nearest'),'r','transparency',.1)
        set(gca,'xscale','log')
%         set(gca,'yscale','log')
        axis([0 150 0 8000])
        end
        if ~isempty(HPC)
%         subplot(3,3,8)
%         imagesc(log(HPC_mat_phase))
%         caxis([0 log(8000)])
%         subplot(3,3,9)
%         imagesc(log(HPC_mat_rate))
%         caxis([0 log(8000)])
        subplot(3,3,8)
%         hold off
        plot(1)
%         boundedline(1:150,fillmissing(HPC_mat_phase(:,wind),'nearest'),fillmissing(HPC_mat_phase_ste(:,wind),'nearest'),'g','transparency',.1)
%         boundedline(1:150,fillmissing(HPC_mat_rate(:,wind),'nearest'),fillmissing(HPC_mat_rate_ste(:,wind),'nearest'),'k','transparency',.1)
        set(gca,'xscale','log')
%         set(gca,'yscale','log')
%         axis([0 150 0 8000])
        end
    end
end 
 

% for i=1:length(d)
%     f = find(LS(:,4)==i);
%     l(i,:)=nanmedian(LS(f,1:3));
%     f = find(HPC(:,4)==i);
%     h(i,:)=nanmedian(HPC(f,1:3));
% %     f = find(LS(:,4)==i);
% %     r(i,:)=nanmedian(REC(f,1:3));
% end
% LS = l;
% HPC=h;
% % REC=r;
% 
% subplot(3,3,7)
% scatter(LS(:,1),LS(:,2),'.b')
% hold on
% scatter(LS(:,1),LS(:,3),'.r')
% axis([0 150 -100 8000])
% 
% subplot(3,3,8)
% scatter(HPC(:,1),HPC(:,2),'.b')
% hold on
% scatter(HPC(:,1),LS(:,3),'.r')
% axis([0 150 -100 8000])
% % 
% % subplot(3,3,9)
% % scatter(REC(:,1),REC(:,2),'.b')
% % hold on
% % scatter(REC(:,1),REC(:,3),'.r')
% % axis([0 150 -100 8000])
% % [x y] = polyfit(depths(f),(lp(f))-(lr(f)),1);
% % y1 = polyval(x,depths(f));
% 
% % LS(isnan(LS(:,2)),:) = [];
% % HPC(isnan(HPC(:,2)),:) = [];
% % REC(isnan(REC(:,2)),:) = [];
% 
% for p = 1:2
%     subplot(3,3,7)
%     [x y] = polyfit(LS(:,1),LS(:,2),p);
%     y1 = polyval(x,1:200);
%     plot(1:200,y1,'b')
%     hold on
%     [x y] = polyfit(LS(:,1),LS(:,3),p);
%     y1 = polyval(x,1:200);
%     plot(1:200,y1,'r')
%     axis([0 150 -100 8000])
% 
%     subplot(3,3,8)
%     [x y] = polyfit(HPC(:,1),HPC(:,2),p);
%     y1 = polyval(x,1:200);
%     plot(1:200,y1,'b')
%     hold on
%     [x y] = polyfit(HPC(:,1),HPC(:,3),p);
%     y1 = polyval(x,1:200);
%     plot(1:200,y1,'r')
%     axis([0 150 -100 8000])
% 
% %     subplot(3,3,9)
% %     [x y] = polyfit(REC(:,1),REC(:,2),p);
% %     y1 = polyval(x,1:200);
% %     plot(1:200,y1,'b')
% %     hold on
% %     [x y] = polyfit(REC(:,1),REC(:,3),p);
% %     y1 = polyval(x,1:200);
% %     plot(1:200,y1,'r')
% %     axis([0 150 -100 8000])
% end
% 
% subplot(3,3,5)
% [x y] = polyfit(LS(:,1),LS(:,2),1);
% y1 = polyval(x,1:200);
% plot(1:200,y1,'.b')
% hold on
% [x y] = polyfit(LS(:,1),LS(:,3),1);
% y1 = polyval(x,1:200);
% plot(1:200,y1,'.r')
% axis([0 150 -100 8000])
% [x y] = polyfit(HPC(:,1),HPC(:,2),1);
% y1 = polyval(x,1:200);
% plot(1:200,y1,'b')
% hold on
% [x y] = polyfit(HPC(:,1),HPC(:,3),1);
% y1 = polyval(x,1:200);
% plot(1:200,y1,'r')
% axis([0 150 -100 8000])
% % [x y] = polyfit(REC(:,1),REC(:,2),1);
% % y1 = polyval(x,1:200);
% % plot(1:200,y1,'.k')
% % hold on
% % [x y] = polyfit(REC(:,1),REC(:,3),1);
% % y1 = polyval(x,1:200);
% % plot(1:200,y1,'k')
% % axis([0 150 -100 8000])
% % end

figure(1);
h = histogram2(HPC(:,1),HPC(:,3),'xbinlimits',[0 100], 'ybinlimits',[0 10000],'numbins',[50 10000/50]);
figure(10)
subplot(3,2,1)
imagesc(h.XBinEdges,h.YBinEdges,h.Values')
caxis([0 20])
xlabel('# of cells in ensemble')
ylabel('MSE')
set(gca,'ydir','normal')
hold on
ind = HPC(:,2)==wind & HPC(:,1) > nCellThresh;
[f] = fit(HPC(ind,1),HPC(ind,3),method{2});
plot(f(1:200),'w')
[a b] = min(abs(f(1:1000)-100));
title([num2str(b) ' cells to reach <15cm error'])
caxis([0 20])
figure(1)
h = histogram2(HPC(:,1),HPC(:,4),'xbinlimits',[0 100], 'ybinlimits',[0 10000],'numbins',[50 10000/50]);
figure(10)
subplot(3,2,2)
imagesc(h.XBinEdges,h.YBinEdges,h.Values')
set(gca,'ydir','normal')
xlabel('# of cells in ensemble')
ylabel('MSE')
hold on
ind = HPC(:,2)==wind & HPC(:,1) > nCellThresh;
[f] = fit(HPC(ind,1),HPC(ind,4),method{2});
plot(f(1:200),'w')
[a bb] = min(abs(f(1:1000)-100));
title([num2str(bb) ' cells to reach <15cm error'])
figure(1)
h = histogram2(LS(:,1),LS(:,3),'xbinlimits',[0 100], 'ybinlimits',[0 10000],'numbins',[50 10000/50]);
figure(10)
subplot(3,2,3)
imagesc(h.XBinEdges,h.YBinEdges,h.Values')
caxis([0 20])
set(gca,'ydir','normal')
xlabel('# of cells in ensemble')
ylabel('MSE')
hold on
ind = LS(:,2)==wind & LS(:,1) > nCellThresh;
[f] = fit(LS(ind,1),LS(ind,3),method{2});
plot(f(1:200),'w')
[a bbb] = min(abs(f(1:1000)-100));
title([num2str(bbb) ' cells to reach <15cm error'])
figure(1)
h = histogram2(LS(:,1),LS(:,4),'xbinlimits',[0 100], 'ybinlimits',[0 10000],'numbins',[50 10000/50]);
figure(10)
subplot(3,2,4)
imagesc(h.XBinEdges,h.YBinEdges,h.Values')
caxis([0 20])
set(gca,'ydir','normal')
xlabel('# of cells in ensemble')
ylabel('MSE')
hold on
ind = LS(:,2)==wind & LS(:,1) > nCellThresh;
[f] = fit(LS(ind,1),LS(ind,4),method{2});
plot(f(1:200),'w')
[a bbbb] = min(abs(f(1:1000)-100));
title([num2str(bbbb) ' cells to reach <15cm error'])

subplot(3,2,5)
bar([b bb bbb bbbb])
ylabel('# of neurons')
xlabel('HPC phase    HPC rate    LS phase     LS rate')

subplot(3,2,6)
caxis([0 20])
colorbar
ylabel('count')

figure(1)
















