%%
clear
close all
load('combnet2.mat')
net = combnet2;
%%
% Enter video file name here with extension : '***.avi'
video='C:\Users\gxs37\Documents\Research Theoretical Biophysics\RBC adhesion\Image processing\Sample-84-3_0_ramp.avi';
vid=VideoReader(video);

%% INITIALIZING FRAME

T=vid.Duration;

% number of frames to be counted in the video
% frames are calculated at equal separation of time
% METHOD 1 - KNOW # of data points
% t=0;
% to=0;
% c=200;
% dt=(T-to)/c;
% numbercells = zeros (1, c);
% time = zeros (1, c);

% Method 2 - KNOW time interval between data points
 t=384.001;
 dt =6 ;

%
vid.CurrentTime = t;
n=1;

frames1=readFrame(vid);

%% IDENTIFYING CELLS Classes

%bigImg = imread('3.png');

bigImg = frames1;
%frames1 = bigImg;

def =0;
nondef = 0;

C = semanticseg(bigImg,net);
B = labeloverlay(bigImg,C);

figure, imshow(B)

D = ( C=="deformable" ) + (C=='nondeformable');
E = bwareaopen(D, 60);
CC = bwconncomp(E);
for i=1:CC.NumObjects
   fin1= mode(C(CC.PixelIdxList{i}));
   fin2 = cellstr(fin1);
   label_str(i) = fin2;
   if strcmp(fin2{1},'deformable')
       label_num(i) = 1;
       colorMap(i) = {'yellow'};
       def = def + 1;
       DefIdArray(def) = i; 
   elseif strcmp(fin2{1},'nondeformable')
       label_num(i) = 2;
       colorMap(i) = {'cyan'};
       nondef = nondef + 1;
       NondefIdArray(nondef) = i;
   end      
end

s = regionprops(CC, 'centroid');
%s = regionprops(E,'centroid');
centroids = cat(1, s.Centroid);
radius=12*ones(size(centroids,1),1);

RGB = insertObjectAnnotation(bigImg,'circle',[centroids(:,1) centroids(:,2) radius],label_str,'FontSize',20,'LineWidth',3,'Color',colorMap,'TextColor','black');
figure
imshow(RGB)
title(['Deformable: ' num2str(def) ', Non-Deformable: ' num2str(nondef) ', Total: ' num2str(nondef + def)] )


%%
for cellid = 1:(def+nondef)
    meanArray(cellid) = mean(bigImg(CC.PixelIdxList{cellid}));
    stdId(cellid) = std2(bigImg(CC.PixelIdxList{cellid}));
end

backgroundPX = mode(bigImg(:));
backgroundDev = std2(bigImg);
CellIdArray = 1: (def+nondef);

%% VIDEO ANALYSIS
vid.CurrentTime = t+dt;
n=n+1;
tracker = input('Enter a number: ');
switch tracker
    case 0
        disp('all cells')
        while hasFrame(vid)
            frames2=readFrame(vid);
            inThisStep =0;
            clear goneCells
            numbercells(1) = length(CellIdArray);
            for cellid = 1:length(CellIdArray)
                avgCellIntensity = mean(frames2(CC.PixelIdxList{CellIdArray(cellid)}));
                if avgCellIntensity < (backgroundPX + 1.5*backgroundDev)
                    inThisStep = inThisStep + 1;
                    goneCells(inThisStep) =(cellid);
                end
            end
            if (exist('goneCells','var'))
                CellIdArray(goneCells) = []; %delete the ones gone
            end
            numbercells(n) = length(CellIdArray);
            %if (numbercells(n)>numbercells(n-1) ) % constraint cell counts to not increase with time
            %    numbercells(n) = numbercells(n-1);
            %end
            %radius of marker circle
            centroids2 = centroids(CellIdArray,:);
            radius2=25*ones(size(centroids2,1),1);
            foundI = insertShape(frames2,'circle',[centroids2(:,1) centroids2(:,2) radius2],'Color','yellow');
            imshow(foundI)
            
            time(n)=t;
            t=t+dt;
            if t<T
                vid.CurrentTime= min(t,T);
            else
                break;
            end
            n=n+1;
        end
        
   
    case 1
        disp('deformable')
        while hasFrame(vid)
            frames2=readFrame(vid);
            inThisStep =0;
            clear goneCells
            numbercells(1) = length(DefIdArray);
            for cellid = 1:length(DefIdArray)
                avgCellIntensity = mean(frames2(CC.PixelIdxList{DefIdArray(cellid)}));
                if avgCellIntensity < (backgroundPX + backgroundDev)
                    inThisStep = inThisStep + 1;
                    goneCells(inThisStep) =(cellid);
                end
            end
            if (exist('goneCells','var'))
                DefIdArray(goneCells) = []; %delete the ones gone
            end
            numbercells(n) = length(DefIdArray);
            %if (numbercells(n)>numbercells(n-1) ) % constraint cell counts to not increase with time
            %    numbercells(n) = numbercells(n-1);
            %end
            %radius of marker circle
            centroids2 = centroids(DefIdArray,:);
            radius2=25*ones(size(centroids2,1),1);
            foundI = insertShape(frames2,'circle',[centroids2(:,1) centroids2(:,2) radius2],'Color','yellow');
            imshow(foundI)
            
            time(n)=t;
            t=t+dt;
            if t<T
                vid.CurrentTime= min(t,T);
            else
                break;
            end
            n=n+1;
        end
        
    case 2
        disp('non-deformable')
        while hasFrame(vid)
            frames2=readFrame(vid);
            inThisStep =0;
            clear goneCells
            numbercells(1) = length(NondefIdArray);
            for cellid = 1:length(NondefIdArray)
                avgCellIntensity = mean(frames2(CC.PixelIdxList{NondefIdArray(cellid)}));
                if avgCellIntensity < (backgroundPX + backgroundDev)
                    inThisStep = inThisStep + 1;
                    goneCells(inThisStep) =(cellid);
                end
            end
            if (exist('goneCells','var'))
                NondefIdArray(goneCells) = []; %delete the ones gone
            end
            numbercells(n) = length(NondefIdArray);
            %if (numbercells(n)>numbercells(n-1) ) % constraint cell counts to not increase with time
            %    numbercells(n) = numbercells(n-1);
            %end
            %radius of marker circle
            centroids2 = centroids(NondefIdArray,:);
            radius2=25*ones(size(centroids2,1),1);
            foundI = insertShape(frames2,'circle',[centroids2(:,1) centroids2(:,2) radius2],'Color','yellow');
            imshow(foundI)
            
            time(n)=t;
            t=t+dt;
            if t<T
                vid.CurrentTime= min(t,T);
            else
                break;
            end
            n=n+1;
        end        
    otherwise
        disp('other value')
end
%

%% Final figure number v time
    figure
    plot (time,numbercells,'r-')
    title('Patient ID:xxxx FlowRate:yyyy')
    xlabel('t (s)')
    ylabel('N (Cell Count)')
    hold off
    A = [time ; numbercells];
    writematrix(A,'Fix_Sample-84-3_0_ramp_ALL.csv')
    
    