



spikes = bz_GetSpikes('getwaveforms',false);
lfp = bz_GetLFP(6);
[b a] = butter(3,[4/625 12/625],'bandpass');
phases = angle(hilbert(FiltFiltM(b,a,double(lfp.data))));

for i=1:length(spikes.times)
for k=1:length(spikes.times{i})
b = round(1250*spikes.times{i}(k));
phasetimes{i}(k) = b;
phaseangles{i}(k) = phases(b);
end
i
end
t = dir('*run*ts*');
ts = load(t.name);
ts(ts==-1)=nan;
pos(:,5) = ts(:,2);
pos(:,2) = ts(:,3);
pos(:,3) = ts(:,2);
pos(:,4) = ts(:,3);
pos(:,1) = ts(:,1) - 93186820;



vel =  [(abs(diff(pos(:,2)))+abs(diff(pos(:,3)))); 0];

amyg = find(spikes.shankID>=5); % rough estimate...

for k=amyg
        p = [phasetimes{k};phaseangles{k}]';
          
        subplot(2,2,3)
        scatter(pos(:,2),pos(:,3),'.k')
        subplot(2,2,4)
        scatter(pos(:,3),pos(:,2),'.k')
        
%         pp = Restrict(spikes.times{k},r.safelaps.all(t,:));
%         for l  =1:length(pp)
%             b = round(1250*pp(l));
%             phasetimes{k}(l) = b;
%             phaseangles{k}(l) = phases(b); 
%         end
  
%         indices = FindInInterval(p(:,1),[pos(1,1) pos(end,1)]);
        if ~isempty(indices)
        for i= 1:length(p) %indices(1):indices(2)
        [a b] = min(abs(p(i,1)-pos(:,1)));
        if vel(b) > nanmean(vel)
        subplot(2,2,1);scatter(pos(b,1),p(i,2),'.k');hold on
%         scatter(pos(b,2),p(i,2)+2*pi,'.k')
%          axis([0 600 -pi pi*3])
        subplot(2,2,2);scatter(pos(b,1),p(i,2),'.k');hold on
%         scatter(pos(b,3),p(i,2)+2*pi,'.k')
%          axis([0 600 -pi pi*3])
        pause(.01)
            end
        end
pause
clf
    end

end