clear 
d = dir('*201*');
clf
LS=[];HPC=[];
LS_mat_phase = nan(150,200);
HPC_mat_phase = nan(150,200);
LS_mat_rate = nan(150,200);
HPC_mat_rate = nan(150,200);
% converge faster with random sampling
% ord = randperm(length(d));
wind = 5; % smoothing window to plot
nCellThresh = 1; % min num cells to be counted

for i=length(d):-1:1
%     cd(d(i).name)
    load(d((i)).name)
for tau = 5
    hpc=[];ls=[];
    if strcmp(positionDecodingMaxCorr_binned_box_mean.region{1},'ls')
        ls = load(d(i).name);
        conditions = unique(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.condition);
    elseif strcmp(positionDecodingMaxCorr_binned_box_mean.region{1},'hpc')
        hpc = load(d(i).name);
        conditions = unique(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.condition);
    end

    for cond = conditions'
%         if length(behavior.events.trialConditions==cond)>=10
        if ~isempty(ls)
        rows = find(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.tau == tau);
        cols = find(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.condition == cond);
        ind = intersect(rows,cols);
        ls_mse_phase(1) = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase(ind));
        ls_mse_phase(2) = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_all(ind));
        ls_mse_phase(3) = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_cos(ind));
        ls_mse_phase = min(ls_mse_phase);
        ls_mse_rate = nanmedian(ls.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_rate(ind));
        ls_nCells = length(ls.positionDecodingMaxCorr_binned_box_mean.region);
        elseif ~isempty(hpc)
        rows = find(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.tau == tau);
        cols = find(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.condition == cond);
        ind = intersect(rows,cols);
        hpc_mse_phase(1) = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase(ind));
        hpc_mse_phase(2) = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_all(ind));
        hpc_mse_phase(3) = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_phase_cos(ind));
        hpc_mse_phase = min(hpc_mse_phase); 
        hpc_mse_rate = nanmedian(hpc.positionDecodingMaxCorr_binned_box_mean.results{end}.mse_rate(ind));
        hpc_nCells = length(hpc.positionDecodingMaxCorr_binned_box_mean.region);
        end

        if ~isempty(ls) & tau == wind
                subplot(3,3,1)
                scatter(ls_nCells,ls_mse_phase,'.b')
                hold on
                scatter(ls_nCells,ls_mse_rate,'.r')
                axis([0 150 -100 8000])
                title('ls cells')
        end
        if ~isempty(hpc) & tau == wind  
                subplot(3,3,2)
                scatter(hpc_nCells,hpc_mse_phase,'.b')
                hold on
                scatter(hpc_nCells,hpc_mse_rate,'.r')
                axis([0 150 -100 11000])
                title('hpc cells')
        end
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
    if ~mod(i,50)
        if ~isempty(LS)
        LS(isnan(LS(:,3)),:) = [];
        end
        if ~isempty(HPC)
        HPC(isnan(HPC(:,3)),:) = [];
        end
        % REC(isnan(REC(:,2)),:) = [];
        subplot(3,3,4)
        hold off
        subplot(3,3,5)
        hold off
        method = {'poly1','exp1'};
        for p = 1:2
            if ~isempty(LS)
            ind = LS(:,2)==wind & LS(:,1) > nCellThresh;
            subplot(3,3,4)
            [f] = fit(LS(ind,1),LS(ind,3),method{p});
%             y1 = polyval(x,1:200);
            plot(f(1:200),'b')
            hold on
            [f] = fit(LS(ind,1),LS(ind,4),method{p});
%             y1 = polyval(x,1:200);
            plot(f(1:200),'r')
            axis([0 150 -100 11000])
            end
            if ~isempty(HPC)
            ind = HPC(:,2)==wind & HPC(:,1) > nCellThresh;
            subplot(3,3,5)
            [f] = fit(HPC(ind,1),HPC(ind,3),method{p});
%             y1 = polyval(x,1:200);
            plot(f(1:200),'b')
            hold on
            [f] = fit(HPC(ind,1),HPC(ind,4),method{p});
%             y1 = polyval(x,1:200);
            plot(f(1:200),'r')
            axis([0 150 -100 11000])
            end
        end
        for i=1:150 % cell num
            if ~isempty(LS)
            f = find(LS(:,1)==i);
            end
            if ~isempty(HPC)
                g = find(HPC(:,1)==i);
            end
            for j=1:200 % tau
                % set minimum number of models to be included?
                if ~isempty(LS)
                ff = find(LS(:,2)==j);
                LS_mat_phase(i,j) = nanmean(LS(intersect(f,ff),3));
                LS_mat_rate(i,j) = nanmean(LS(intersect(f,ff),3));
                end
                if ~isempty(HPC)
                gg = find(HPC(:,2)==j);
                HPC_mat_phase(i,j) = nanmean(HPC(intersect(g,gg),3));
                HPC_mat_rate(i,j) = nanmean(HPC(intersect(g,gg),3));
                end
            end
        end
        subplot(3,3,6)
        imagesc(log(LS_mat_phase))
        ylabel('cell #')
        xlabel('smoothing window')
        subplot(3,3,7)
        imagesc(log(LS_mat_rate))
        subplot(3,3,8)
        imagesc(log(HPC_mat_phase))
        subplot(3,3,9)
        imagesc(log(HPC_mat_rate))
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