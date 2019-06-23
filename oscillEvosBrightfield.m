%% SET FOLDER AND FILENAMES, SET PROCESSING PARAMETERS, INITIALIZE VARIABLES
clear; close all; clc;

date = input('Please enter the date the images were recorded (the folder name): ','s');
beacon = input('Please enter the beacon you want to analyze: ','s');
base = 'file path here'; %MANUALLY CHANGE DIRECTORY
imdir= [base date '\Beacon-' beacon '\'];
cd(imdir); %so script and images don't have to be in same directory
images= dir([imdir '\*.tif']); %MANUALLY CHANGE FILE TYPE


n = length(images); %number of images in folder 
adjust = input('Please enter the number of images in the folder that are not to be processed (e.g. scale and reference): ');

step = 1; %MANUALLY CHANGE depending on how often want to analyze a frame - 1 means every frame
j = ceil((n-adjust)/step); %this is number of frames that will be processed

% I is list of filenames in numerical order
I{n-adjust,1} = [];
imname{n-adjust,1} = [];

for i = 1:n-adjust;
    if i<10
        imname{i} = ['Scene1Interval00' num2str(i)]; %MANUALLY CHANGE BASE NAME %only need i-1 if images numbered starting at 0
        I{i} = [imdir imname{i} '.tif'];  %MANUALLY CHANGE FILE TYPE
        
    elseif i<100
        imname{i} = ['Scene1Interval0' num2str(i)]; %MANUALLY CHANGE BASE NAME %only need i-1 if images numbered starting at 0
        I{i} = [imdir imname{i} '.tif'];  %MANUALLY CHANGE FILE TYPE
        
    else
        imname{i} = ['Scene1Interval' num2str(i)]; %MANUALLY CHANGE BASE NAME %only need i-1 if images numbered starting at 0
        I{i} = [imdir imname{i} '.tif'];  %MANUALLY CHANGE FILE TYPE
        
    end
        
end

imnum = 1;
info = imfinfo(images(imnum).name);

%Initialize cell array for saving values
%Note: use cell arrays instead of basic array because can have multiple data types
M{j,9}=[]; %MANUALLY CHANGE SECOND NUMBER DEPENDING ON VALUES SAVING
%% USE REFERENCE IMAGE AND FIRST IMAGE TO ESTABLISH BACKGROUND SUBTRACTION
% good = 0;
% msg = 'Thats not good. Make sure your reference image is correct and try again.';
% while (good==0)
%     imdir= [base date 'ref\well' well];
%     ref = imread([imdir '\ref_0.tif']); %MANUALLY CHANGE if reference named something else
%     first = imread(I{1});
%     figure, subplot(2,2,1), imshow(ref), subplot(2,2,2), imshow(first);
%     sub = first - ref;
%     subplot(2,2,3),imshow(sub);
% 
%     choice = questdlg('Is the background subtraction satisfactory?','Settings','Yes','No', 'Yes'); %ask if image needs to be adjusted
%     % handle response
%     switch choice
%         case 'Yes'
%             good = 1;
%         case 'No'
%             error(msg) %terminate script if 'no' because no mechanism for modifying background subtraction
%             
%     end
% end
% close all;
%% ESTABLISH PROCESSING SETTINGS USING FIRST IMAGE
loop = 0;
while (loop==0)
    im = imread(I{1}); % import the first image and get size
    [height,width] = size(im);
    
    %im = im - ref; %background subtracted
    figure(1); imshow(im) %display first image
    
    %create settings box
    prompt = {'enter cropping window size:','enter area filter size:'}; %leave crop for flexibility to get better threshold even though only one well per image
    dlg_title = 'settings';
    num_lines = 1;
    def = {'300','100'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    % answers to settings box
    idealcrop = str2double(answer(1)); %cropping window size; needs to be specified manually;
    af = str2double(answer(2)); % area filter size (how small of objects to clear)

    %input specified settings for cropping image around ginput
    [c0,r0] = ginput(1); r0=round(r0); c0=round(c0);% row/column, origin in upper left corner, ginput works on the figure not the image
    wx=min([idealcrop width-c0 c0-1]);
    wy=min([idealcrop height-r0 r0-1]);
    im_crop = im((r0-wy):(r0+wy),(c0-wx):(c0+wx)); %crops and makes new image

    figure(2); subplot(1,2,1); imshow(im_crop); %displays cropped image

    % thresholding specifications for locating object
    level = graythresh(im_crop); % compute global threshold value
    new_im = im2bw(im_crop,level); %threshold cropped image using calculated threshold
    imc = imcomplement(new_im);
    bw = imclose(imc,ones(5,5));
    bw2 = imfill(bw,'holes'); %fill holes
    bw3 = bwareaopen(bw2,af);
    bw4 = imclearborder(bw3); %clear objects on border -- this creates an error if the tracking fails
    
    figure(2);subplot(1,2,2); imshow(bw4); %displays thresholded cropped image
    
    choice = questdlg('Is the image satisfactory?','Settings','Yes','No', 'Yes'); %ask if image needs to be adjusted
    % handle response
    switch choice
        case 'Yes'
            loop = 1;
        case 'No'
            loop = 0;
    end

end

close all;
%% PROCESS DESIRED FRAMES FOR ENTIRE RUN

clear x y stats
start=1; 
finish=n-adjust;

x=zeros(j,1);  %pre-allocate space for matrices
y=zeros(j,1);
com=zeros(j,2);

count = 1; %keeps track of how many frames have been processed
for i = start:step:finish; 
    im = imread(I{i});
    dimr = size(im, 1); dimc = size(im, 2);
    r   = r0;   c   = c0;
    
    %im = im - ref; %background subtracted
    imc = im;
%   imc = imcomplement(imc); % use if worm darker than background

    imc= double(imc);
    imc=(imc*255/max(max(imc)));
    imc=uint8(imc); %converts image to grayscale(?)

    wx=min([idealcrop width-c0 c0-1]);
    wy=min([idealcrop height-r0 r0-1]);
    im_crop = imc((r-wy):(r+wy),(c-wx):(c+wx)); %crops and makes new image
    

    level = graythresh(im_crop); % compute global threshold value
    new_im = im2bw(im_crop,level); %threshold cropped image using calculated threshold
    imc = imcomplement(new_im);
    bw = imclose(imc,ones(5,5)); %close image to improve segmentation around edge
    bw2 = imfill(bw,'holes'); %fill holes
    bw3 = bwareaopen(bw2,af); %clear small objects (artifacts of thresholding)
    bw4 = imclearborder(bw3); %clear objects on border -- this creates an error if the tracking fails

    cc = bwconncomp(bw4); %find connected components in image to identify objects
    
    stats=regionprops(cc,'Area','Centroid', 'MajorAxisLength','MinorAxisLength', 'Orientation'); 
    %Find area and other relevant properties of objects found.

    % Determine which of the remaining objects has the largest area
    % This is the relevant object (the aggregate)
    temp = 0;
    for j = 1:length(stats)
        if stats(j).Area > temp 
            temp = stats(j).Area;
            label = j;
        end
    end

    %Plot the thresholded object onto the original image for comparison.
%         [B,L] = bwboundaries(bim6,'noholes');
%         figure(3);imshow(im_crop);
%         hold on    
%         
%         for k = 1:length(B);
%             boundary = B{k};
%             plot(boundary(:,2), boundary(:,1), 'm', 'LineWidth', 2)
%         end
%         
%         hold off;

    com(count,:) = stats(label).Centroid; %for tracking object
    A=stats(label).Area; %tells you the area in pixel^2
    MMin=stats(label).MinorAxisLength; %in pixel
    MMax=stats(label).MajorAxisLength; %in pixel
    R = MMin/MMax; %aspect ratio, closer to 1 is more circular
    Angle = stats(label).Orientation; %in degrees, -90 to 90
    Centroid = stats(label).Centroid; %for drawing ellipse


    %Draw fitted ellipse and major/minor axes periodically
    if i == start || i == finish || mod(count,72) == 0 %%MANUALLY CHANGE depending on desired frequency
        t = linspace(0,2*pi,50);
        a = MMax/2;
        b = MMin/2;
        Xc = Centroid(1);
        Yc = Centroid(2);
        phi = deg2rad(-Angle);
        ellipx = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        ellipy = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);

        figure, imshow(im_crop)
        hold on;
        plot(ellipx, ellipy,'m', 'LineWidth', 2);

        % Draw lines along the axes
        xMajor = Xc + [-1 1]* a*cos(phi);
        yMajor = Yc + [-1 1]* a*sin(phi);
        line(xMajor,yMajor,'Color','m','LineWidth', 2);

        xMinor = Xc + [-1 1]* b*cos(phi - pi/2);
        yMinor = Yc + [-1 1]* b*sin(phi - pi/2);
        line(xMinor,yMinor,'Color','m','LineWidth', 2);

        hold off;
    end

    All = [A, MMin, MMax, R, Angle, Centroid(1), Centroid(2)];
    Bll = num2cell(All);
    M(count,:) = [count I{i} Bll]; %store measured values for current image

    xcrop = com(count,1);
    ycrop = com(count,2);
    x(count)= c - wx - 1 + xcrop; % position in current original image
    y(count)= r - wy - 1 + ycrop; % position in current original image
    r0 = round(y(count)); % row, for positioning the next cropping window
    c0 = round(x(count)); % column
    
    count = count + 1;
end
%% SAVE RESULTS AS MAT FILE
saveBase = '\\planaria.ucsd.edu\HYDRA\Kate Khazoyan\Evos\'; %MANUALLY CHANGE SAVE DIRECTORY
savedir= [saveBase date '\'];
cd(savedir);
fileBase = 'oscill-';
filename= [fileBase beacon];
save(filename,'M'); 
%% PLOT RESULTS AND SAVE SIMPLIFIED RESULTS MATRIX
%Get input for plotting parameters and make time vector for x axis
hoursImaged = input('Please enter the number of hours of imaging data: ');
conversion = input('Please enter the pixels/mm length conversion for these images: ');
areaConversion = conversion^2;
time = linspace(0,hoursImaged,count-1); % convert time axis from frames to time

%Turn measurements in M into plottable vectors
areaPix = cell2mat(M(1:count-1,3));
radiusPix = sqrt(areaPix./pi); 
areaMM = areaPix./areaConversion; %convert area from pixels^2 to mm^2
radiusMM = sqrt(areaMM./pi); %convert area to average radius in mm -- how done in Soriano et al 2009
radiusUM = radiusMM.*1000; %convert radius to um
ratio = cell2mat(M(1:count-1,6));
angle = cell2mat(M(1:count-1,7));

%Plot area/radius, aspect ratio, and angle vs. time, and save plots
%figure,plot(time,areaMM, '.:'), title(['Area vs. Time ' date ' ' well]), xlabel('Time (hours)'), ylabel('Area(mm^2)')
figure,plot(time,radiusUM, '.:'), title(['Radius vs. Time ' date ' ' beacon]), xlabel('Time (hours)'), ylabel('Radius(um)')
% % figure,plot(time,radiusPix, '.:'), title(['Radius vs. Time ' date ' ' well]), xlabel('Time (hours)'), ylabel('Radius(pixels)')
savefig([beacon '-radius.fig']);
figure,plot(time,ratio,'.:'), title(['Aspect Ratio vs. Time ' date ' ' beacon]), xlabel('Time (hours)'), ylabel('Aspect Ratio(1 = circle)')
savefig([beacon '-ratio.fig']);
figure,plot(time,angle,'.:'), title(['Orientation Angle vs. Time ' date ' ' beacon]), xlabel('Time (hours)'), ylabel('Angle(-90 to 90 degrees)')
savefig([beacon '-angle.fig']);

M_mod=zeros(count-1,4);
M_mod(:,1)= time;
M_mod(:,2)= radiusUM;
M_mod(:,3)= ratio;
M_mod(:,4)= angle;

fileBase = 'oscill-mod-';
filename= [fileBase beacon];
save(filename,'M_mod'); 