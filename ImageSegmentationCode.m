function iNetoCounts2
%load('combnet2.mat')
%net = combnet2;
%%
%load('combnet2.mat')
%net = combnet2;

bigImg = imread('33 Normoxia Count 1006.jpg');
figure, imshow(bigImg)

%%

def =0;
nondef = 0;

%%
tic 
fun = @(block_struct) segNbinarize(block_struct.data);
doubimg = blockproc(bigImg,[224 224],fun);
toc
%%
C = categorical(doubimg);
%C = semanticseg(bigImg,net);
B = labeloverlay(bigImg,C);

%%
imshow(B)

%%
D = ( doubimg==1 ) + (doubimg==2);
%D = ( C=="deformable" ) + (C=='nondeformable');
D2 = bwareaopen(D, 60);
E = D2 - bwareafilt(D2,[400 5000]);
E=logical(E);
CC = bwconncomp(E);
for i=1:CC.NumObjects
   fin1= mode(C(CC.PixelIdxList{i}));
   fin2 = cellstr(fin1);
   label_str(i) = fin2;
   if strcmp(fin2{1},'1')
%    if fin2{1}==1
       label_num(i) = 1;
       colorMap(i) = {'yellow'};
       def = def + 1;
   elseif strcmp(fin2{1},'2')
%    elseif fin2{1}==2
       label_num(i) = 2;
       colorMap(i) = {'cyan'};
       nondef = nondef + 1;
   end      
end

%%

%Sxy = regionprops(CC,'Centroid')
s = regionprops(E,'centroid');
centroids = cat(1, s.Centroid);
%%
radius=12*ones(size(centroids,1),1);

RGB = insertObjectAnnotation(bigImg,'circle',[centroids(:,1) centroids(:,2) radius],label_str,'FontSize',20,'LineWidth',3,'Color',colorMap,'TextColor','black');
figure
imshow(RGB)
title(['Deformable: ' num2str(def) ', Non-Deformable: ' num2str(nondef) ', Total: ' num2str(nondef + def)] )

%%
end
function [bwimg] = segNbinarize(img)
if ~exist('transferunet1','var')
    load('C:\Users\gxs37\Documents\Research Theoretical Biophysics\RBC adhesion\Image processing\Deep Learning\transferunet1.mat')
    net = transferunet1;
end

%net = combnet2;
segimg = semanticseg(img,net);
bwimg = zeros(size(segimg));
bwimg(segimg=="deformable")=1;
bwimg(segimg=="nondeformable")=2;
end
