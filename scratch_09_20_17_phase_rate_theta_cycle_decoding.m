


for pre = 1:39
    for cond = 1:10
        [pks locs] = findpeaks(phasetrains{cond});
        response = spktrains{cond}(pre,:);
        spkr = find(response==1);
        for cell = 1:74
        p = peers{cond}(cell,:);
        spk = find(p==1);
        if ~isempty(spk)
        count=1; clear num pha lsnum lspha
        for i=1:length(locs)-1
        c = Restrict(spk,[locs(i) locs(i+1)]);
        cr = Restrict(spkr,[locs(i) locs(i+1)]);
        if ~isempty(c)
        num(count) = length(c);
        pha(count) = circ_mean(phasetrains{cond}(c)');
        lsnum(count) = length(Restrict(spkr,[locs(i) locs(i+1)]));
        if ~isempty(cr)
            lspha(count) = circ_mean(phasetrains{cond}(cr)');
        else
            lspha(count) = 0;
        end
        count=1+count;
        end
        end
        
% u = unique(num);
% for i=1:length(u)
% f =find(u(i) == round(unique(round(num./2))));
% rrs(i) = std(lsnum(f));
% rr(i) = mean(lsnum(f));
% end
% subplot(2,2,1)
% boundedline(u,rr,rrs)
% u = unique(round(-pi:.5:pi));
% for i=1:7
% f =find(u(i) == round(pha));
% rp(i) = mean(lsnum(f));
% rs(i) = std(lsnum(f));
% end
% subplot(2,2,2);
% boundedline(u,rp,rs)
% pause(.01)
% clf
% clear r rr rs rrs

        if count > 5
        [phaseLS(pre,cell,cond) b] = circ_corrcl(pha,lsnum);
        [countLS(pre,cell,cond) b] = corr(num',lsnum');
        for iter = 1:10
        r= randperm(length(pha));
        train = r(1:round(length(r)/2));
        test = r(round(length(r)/2)+1:end);
        [res dev] = glmfit(cos(pha(train)),lsnum(train),'normal');
        [resn devn] = glmfit(num(train),lsnum(train),'normal');
        yfit = glmval(res,pha(test),'identity');
        yfitn = glmval(resn,num(test),'identity');
        mse(pre,cond,cell,iter) = mean((yfit-lsnum(test)').^2);
        msen(pre,cond,cell,iter) = mean((yfitn-lsnum(test)').^2);
        
        [res dev] = glmfit([cos(pha(train)); sin(pha(train))]',lspha(train),'normal');
        [resn devn] = glmfit(num(train),lspha(train),'normal');
        
        rando = [cos(pha(train)); sin(pha(train))]'; rando = rando(randperm(length(rando)),:);
        [res_shuffle dev_shuffle] = glmfit(rando,lspha(train),'normal');
        
        yfit = glmval(res,[cos(pha(test)); sin(pha(test))]','identity');
        yfitn = glmval(resn,num(test),'identity');
        yfit_shuffle = glmval(res_shuffle,[cos(pha(test)); sin(pha(test))]','identity');
        
        mse_lsphase(pre,cond,cell,iter) = mean((yfit-lspha(test)').^2);
        msen_lsphase(pre,cond,cell,iter) = mean((yfitn-lspha(test)').^2);
        mse_lsphase_shuffle(pre,cond,cell,iter) = mean((yfit_shuffle-lspha(test)').^2);
        end

        errorbar(mean(squeeze(mse_lsphase(pre,cond,cell,:))),std(squeeze(mse_lsphase(pre,cond,cell,:))),'g')
        hold on
        errorbar(mean(squeeze(msen_lsphase(pre,cond,cell,:))),std(squeeze(msen_lsphase(pre,cond,cell,:))),'r')
        errorbar(mean(squeeze(mse_lsphase_shuffle(pre,cond,cell,:))),std(squeeze(mse_lsphase_shuffle(pre,cond,cell,:))),'k')
        pause(.1)
        hold off
        end
        end
        end
    end
end







