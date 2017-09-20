


for pre = 1:39
    for cond = 1:10
        [pks locs] = findpeaks(phasetrains{cond});
        response = spktrains{cond}(pre,:);
        spkr = find(response==1);
        for cell = 1:74
        p = peers{cond}(cell,:);
        spk = find(p==1);
        if ~isempty(spk)
        count=1; clear num pha lsnum
        for i=1:length(locs)-1
        c = Restrict(spk,[locs(i) locs(i+1)]);
        if ~isempty(c)
        num(count) = length(c);
        pha(count) = circ_mean(phasetrains{cond}(c)');
        lsnum(count) = length(Restrict(spkr,[locs(i) locs(i+1)]));
        count=1+count;
        end
        end
        if count > 1
        [phaseLS(pre,cell,cond) b] = circ_corrcl(pha,lsnum);
        [countLS(pre,cell,cond) b] = corr(num',lsnum');
        for iter = 1:10
        r= randperm(length(pha));
        train = r(1:round(length(r)/2));
        test = r(round(length(r)/2):end);
        [res dev] = glmfit(pha(train),lsnum(train),'normal');
        [resn devn] = glmfit(num(train),lsnum(train),'normal');
        yfit = glmval(results,pha(test),'identity');
        yfitn = glmval(results,num(test),'identity');
        mse(pre,cond,cell,iter) = mean((yfit-lsnum(test)').^2);
        msen(pre,cond,cell,iter) = mean((yfitn-lsnum(test)').^2);
        end
        end
        end
        errorbar(mean(squeeze(mse(pre,cond,cell,:))),std(squeeze(mse(pre,cond,cell,:))),'g')
        hold on
        errorbar(mean(squeeze(msen(pre,cond,cell,:))),std(squeeze(msen(pre,cond,cell,:))),'r')
        pause(.1)
        hold off
        end
    end
end
