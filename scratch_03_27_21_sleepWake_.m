clear all
d = dir('_*mat');
c = 1;
clf
% for thresh = .2:.01:1
thresh = 0.95;
    for i=1:length(d)
        load(d(i).name,'ccc','idx','count','isig_levels','states')
        ccc(end:length(count),:)=nan;
        
%         subplot(4,2,i)
%         plot(states>.95)
%         hold on
%         plot(states<-.95)
        
        id  = states > thresh;
        idd = states < -thresh;
        id = intersect(idx,find(diff(id)==1));
        idd = intersect(idx,find(diff(idd)==1));
        for j=1:length(id)
            if id(j) > 50 & id(j) < length(count)-50
                cor_trans_w{i}(j,:,:) = ccc(id(j)-50:id(j)+50,:);
            end
        end
        for j=1:length(idd)
            if idd(j) > 50 & idd(j) < length(count)-50
                cor_trans{i}(j,:,:) = ccc(idd(j)-50:idd(j)+50,:);
            end
        end
        
        cor_w(i,:)=nanmean(ccc(intersect(idx,find(states>thresh)),:));
        cor(i,:)=nanmean(ccc(intersect(idx,find(states<-thresh)),:));
    end
%     c=c+1;
    
    
    
for i=1:8
    for j=1:101
        a(i,j) = min(nanmean(cor_trans_w{i}(:,j,41:44)));
        b(i,j) = min(nanmean(cor_trans{i}(:,j,41:44)));
    end
end

boundedline(1:101,median(b),std(b),'k')
hold on
boundedline(1:101,median(a),std(a),'r')
title('red/NREM-to-wake; black/wake-to-NREM')
% plot(mean(b),'k')
ylim([-.4 0])
xline(51)
% end

for i = 1:8
   an(i) = min(cor(i,41:45));
   an_w(i) = min(cor_w(i,41:45)); 
end
figure
plot([an;an_w],'k')
hold on
plot([1 2],[nanmean(an),nanmean(an_w)],'r')
xlabel('sleep                                  wake')
ylabel('correlation')
ylim([-.45 0])
