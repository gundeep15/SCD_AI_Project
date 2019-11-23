%% Biotheory Group at Case Western Reserve University

% Line 4 and 6 should be changed, depending on image directory. Change your path to 
% your directroy containing ImageName = '' (e.g., microfluidic channel images). In
% addition, make sure that your neural networks are loaded: SCD_WB_net and
% SCD_ResNet. 
% 
ImageName = '1 Normoxia Count 677' % Image path 
Img = imread(fullfile(strcat(ImageName, '.jpg'))); % Microfluidic Channel corresponding Image
path = 'C:\Users\npral\Google Drive\BiotheoryGroup\SickleCells\SickleCellDetection_Counting\CroppedMicroFluidicImages\'

%% Do not edit any line of code beneath the following dotted line
%% -----------------------------------------------------------------
%%

alpha = 50; % vertical slices
beta = 20; % horizontal slices 
Images = struct


yy = length(Img(:, 1, 1)) / beta; % batch slices' vertical length
xx = length(Img(1, :, 1)) / alpha; % batch slices' horizontal length
z = 0;
zz = 1;

% Line 17 is optional
% figure, 
for ii = 1 : beta
    
   for jj = 1 : alpha
       
       rect = [(0 + (jj - 1) * xx) (0 + (ii - 1) * yy) ...
           xx-1 yy-1];
       
       Images(ii,jj).ChannelBatches = imcrop(Img, rect);
       Images(ii,jj).ChannelBatches = imresize(Images(ii,jj).ChannelBatches,[224, 224]);
       
       z = z + 1;
       zz = zz + 1;
      
      % Optional -- save channel batches  
      
      
      path3 = strcat(path, ImageName, '\Batch_', num2str(z), '.png');
      imwrite(Images(ii, jj).ChannelBatches, path3);
   end 
end
%% Neural Network Mask Prediction + Sickle Cell Extraction:
 clear SegmentVars
SegmentVars = struct;

def = 0;
nondef = 0;

% Line 44 and 62 are optional
% figure, 

for ii = 1 : beta
    for jj = 1 : alpha
        
        % make sure to load "net" into your workspace
        % Mask prediction 
        SegmentVars(ii,jj).Prediction = semanticseg(Images(ii,jj).ChannelBatches, SCD_WB_Net);
        
       
        mask1 = SegmentVars(ii,jj).Prediction == 'deformable'  ;% Deformable Binary Image
        mask2 = SegmentVars(ii,jj).Prediction == 'nondeformable' ;% Nondeformable Binary Image
        mask = mask1 + mask2 ; % Deformable + Nondeformable Binary Image
        mask = bwareaopen(mask, 60); % Removing small binarized pixels 
        
        % Illustration
        E = labeloverlay(Images(ii,jj).ChannelBatches, ...
             mask);
         
%       Line 61 is optional
      imshow(E) 
        
        % Categoricalize Predicted Masks
        SegmentVars(ii, jj).categorical = categorical(SegmentVars(ii,jj).Prediction);
        

    end
end

%%

for ii = 1 : beta 
    for jj = 1 : alpha
        
        ii 
        jj
        SegmentVars(ii,jj).Def_Nondef_Labels = ( SegmentVars(ii,jj).categorical == 'deformable')  + (SegmentVars(ii,jj).categorical == 'nondeformable');
        SegmentVars(ii,jj).BwAreaOpen = bwareaopen(SegmentVars(ii,jj).Def_Nondef_Labels, 60);
        SegmentVars(ii,jj).BWconn = bwconncomp(SegmentVars(ii,jj).BwAreaOpen);
        
        for kk = 1 : SegmentVars(ii,jj).BWconn.NumObjects
            
            
            
            fin1 = mode(SegmentVars(ii,jj).categorical((SegmentVars(ii,jj).BWconn.PixelIdxList{kk})));
            fin2 = cellstr(fin1);
            label_str(ii) = fin2;
            
            % If deformable cell is detected, then add one to 'def'
            if strcmp(fin2{1}, 'deformable')
                
                label_num(ii) = 1;
                colorMap(ii) = {'yellow'};
                def = def + 1;           
             
            % If nondeformable cell is not detect, then add one to 'nondef'   
            elseif strcmp(fin2{1}, 'nondeformable')
                
                label_num(ii) = 2;
                colorMap(ii) = {'cyan'};
                nondef = nondef + 1;
                
            end
        end
        
    end
end
%% Locating Centroids: 

% Creating a matrix based on each cell's centroid coordinates
for ii = 1 : beta
    for jj = 1 : alpha
         
        Centroid(ii,jj).def_nondef = regionprops(SegmentVars(ii,jj).BwAreaOpen, 'centroid');
        Centroid(ii,jj).centroids = cat(1, Centroid(ii, jj).def_nondef);
  
    end
end
%% Sanity Check:

z = 0;

for ii = 1 : beta
    for jj = 1 : alpha
        
        z = z + length(Centroid(ii,jj).centroids);
        
    end
end

disp('Sanity Check: ')
Difference_Cells = z - nondef - def


%% 

temp1 = 0;
xmin = 0;
ymin = 0;

for ii = 1 : beta
    for jj = 1 : alpha
        
        for kk = 1 : length(Centroid(ii,jj).centroids)
            
            temp1 = temp1 + 1;
            coord = struct2cell(Centroid(ii,jj).centroids(kk));
            coord2 = cell2mat(coord);
            xmin = round(coord2(1,1) - 10);
            ymin = round(coord2(1,2) - 15);
            width = 31; 
            height = 31;
            rect = [xmin, ymin, width, height];
            J = imcrop(Images(ii,jj).ChannelBatches, rect); 
            J = imresize(J, [224, 224]);
            
           imwrite(J,strcat(path, ImageName, '\CroppedCells\Frame_', num2str(temp1),'.png'));            ii;

            
        end  
    end 
end

%%

pathNew = strcat(path, ImageName, '\CroppedCells')

zz_def = 0; 
zz_nondef = 0;
zz_UFO = 0; 

pathNew2 = strcat(path,ImageName,'\CorrectedCell\UFO')

for ii = 1 : (def + nondef)
    
    CellPath = strcat(pathNew, '\Frame_', num2str(ii), '.png');
    
    Img = imread(CellPath);
    imshow(Img);
    YPredicted = classify(SCD_ResNet, Img)
    
    if YPredicted == 'Resized_nondeformable'
        
        zz_nondef = zz_nondef + 1;
        imwrite(Img, strcat(path, ...
            ImageName, '\CorrectedCell\nondeformable\NondeformableCell_', num2str(zz_nondef),'.png'));
        
    elseif YPredicted == 'Resized_deformable'
        
        zz_def = zz_def + 1;
        
        imwrite(Img,strcat(path, ... 
            ImageName, '\CorrectedCell\deformable\DeformableCell_', num2str(zz_nondef),'.png')); 

    elseif YPredicted == 'Resized_UFO' 
        
       zz_UFO = zz_UFO + 1;
       
       imwrite(Img,strcat(path, ...
           ImageName, '\CorrectedCell\UFO\UFOCell_', num2str(zz_nondef),'.png'));
            

    end
    
end
%%
ImageName
disp('Def = '); disp(zz_def)
disp('nondef = '); disp(zz_nondef)
disp('UFO = '); disp(zz_UFO)
disp('Total = '); disp(zz_def + zz_nondef)

