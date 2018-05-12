% so the main things youll probably want to do are smooth with imgaussfilt, 
%threshold, remove extra blobs with bwareaopen, 
%then find the 1D skeleton with bwmorph(img,'skel',inf)



tifs = dir('*tif');
image = imread(tifs(1).name);
[height, width] = size(image);

% Initialize grid for locating centroid of LED
[X,Y] = meshgrid(1:width,1:height);


imshow(image);
[x,y] = ginput;
inArea = inpolygon(X(:),Y(:),x,y);
inArea = double(reshape(inArea,[height width]));
image(~inArea) = 0;

clf,imshow(image)

image2 = imgaussfilt(image,4);
image2 = image2>40;
image3 = bwmorph(image2,'skel',inf);
clf, imshow(image3)

image = image(:,4500:end);
bw = image~= 0;
for i =1:10
    image2 = imgaussfilt(image,5*i);
    image2 = image2>30;
    image3 = bwmorph(image2,'skel',inf);
    subplot(5,3,i);
    image4 = bwskel(image3,'minbranchLength',1000);
    imshow(image3)
%     image2(~bw) = Inf;
%     L = watershed(image2);
%     L(~bw) = 0;
%     rgb = label2rgb(L,'jet',[.5 .5 .5]);
%     figure(1)
%     subplot(5,3,i)
%     imshow(rgb)
%     title('Watershed transform of D')
    pause(.1)
end

