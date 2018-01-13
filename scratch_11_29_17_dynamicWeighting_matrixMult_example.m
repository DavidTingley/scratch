sin(0:.01:pi);
sin(0:.01:pi*100);
d = sin(0:.01:pi*100);
for i =1:50
for i =1:315
mat(i,:)= circshift(d,i);
end
whos mat
imagesc(mat)
for i =1:315
mat2(i,:) = ones(31416,1);
for k=1:31416
mat2(i,k) = mat2(i,k) + i*k;
end
end
whos mat2
imagesc(mat2)
imagesc(mat2.*mat)
whos mat2
imagesc(mat2)
max(mat2(:))
for i =1:315
mat2(i,:) = ones(31416,1);
for k=1:31416
mat2(i,k) = mat2(i,k) + i*k-(9896041-i*k);
end
end
imagesc(mat2)
caxis
for i =1:315
mat2(i,:) = ones(31416,1);
for k=1:31416
mat2(i,k) = mat2(i,k) + i*k+(9896041-i*k);
end
end
imagesc(mat2)
caxis
dd = [1:31416/2 31416/2:-1:1];
plot(d)
plot(dd)
for i =1:315
mat2(i,:) = circshift(dd,31416/315*i);
end
for i =1:315
mat2(i,:) = circshift(dd,round(31416/315*i));
end
imagesc(mat2)
dd = [31416:-1:1];
for i =1:315
mat2(i,:) = circshift(dd,round(31416/315*i));
end
imagesc(mat2)
dd = [1:31416/2 31416/2:-1:1];
whos dd
for i =1:315
mat2(i,:) = circshift(dd,31416/2+round(31416/315*i));
end
imagesc(mat2)
imagesc(mat2.*mat)
imagesc(mat2)
imagesc(mat2.*mat)
dd = [1:31416 31416:-1:1];
for i =1:315
mat2(i,:) = circshift(dd,31416/2+round(31416/315*i));
end
clear mat2
for i =1:315
mat2(i,:) = circshift(dd,31416/2+round(31416/315*i));
end
whos mat2
mat2 = mat2(:,1:31416);
imagesc(mat2)
imagesc(mat2.*mat)
dd = [1:31416/2 31416/2:-1:1];
for i =1:315
mat2(i,:) = circshift(dd,31416/2+round(31416/315*i));
end
imagesc(mat2)
imagesc(mat2.*mat)
for i=1:31416
m = mat2.*mat;
for i=1:31416
[a b] = max(m(:,i));
ind(i) = b;
end
plot(ind)
plot(ind,'.')
imagesc(mat2.*mat)
whos d
plot(d)
plot(d,ind)
plot(d,ind,'.')
plot(d,ind)
plot(d,ind,'.')
plot(d,d,'.')
plot(ind,d,'.')
for i=1:31416
[a b] = max(m(:,i));
ind(i) = b;
w(i) = a;
end
plot(ind,d,'.')
plot(w,d,'.')
plot(w,ind,'.')
scatter3(w,ind,d,'.')
min(d)
max(d)
clf
imagesc(mat2.*mat)
plot(w,ind,'.')
plot(ind,d,'.')
whos
plot(ind)
plot(ind,'.')
histogram(w)
plot(ind(w>10000),'.')
plot(-ind(w>10000),'.')
plot(d,-ind(w>10000),'.')
plot(d(w>10000),-ind(w>10000),'.')
plot(d(w>10000),'.')
hold on
plot(d(w>10000)+1,'.')
clf
plot(d(w>10000),'.')
hold on
plot(d(w>10000)+2,'.')
plot(d(w>10000)+2,'.k')
plot(d(w>10000),'.k')
plot(d(w>20000),'.r')
plot(d(w>15000),'.r')
clf
plot(d(w>15000),'.r')
hold on
plot(d(w>15000)+2,'.r')
subplot(2,2,1
subplot(2,2,1)
imagesc(mat2.*mat)
imagesc(mat2)
subplot(2,2,2)
imagesc(mat)
subplot(2,2,3)
imagesc(mat2.*mat)
subplot(2,2,3)
subplot(2,2,4)
plot(d(w>15000)+2,'.r')
hold on
plot(d(w>15000),'.r')
xlabel('time')
ylabel('ca3                ca1')
ylabel('ca3                                 ca1')
xlabel('time/space')
ylabel('ca1 theta phase')
max(w)
mean(w)
std(w)
mean(w)+3*std(w)
mean(w)+2*std(w)
figure
imagesc(mat2)ex
for i =1:315
mat2(i,:) = [1:31416]./(316-i);
end
