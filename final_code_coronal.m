%% Medical Imaging
% This script visualizes 3 planes of DICOM slices (axial, coronal, sagittal)
% with a focus on the coronal slices. The workflow includes:
% Step 1: Automatic or manual cropping of a coronal slice.
% Step 2: Interactive selection of pulmonary branches to be removed.(1 clic + enter)
% Step 3: Interactive selection of the two lungs.(2 clicks + enter)
% Step 4: Loop to segment lungs in each slice
% Step 5: Area calculation
% Step 6: Volume viewer


%% Close figures, clear workspace, and command window
close all;
clear;
clc;

%% Check if the function file exists
if exist('readDCMfolder.m', 'file') ~= 2
    error('The function readDCMfolder.m is not found in the current directory.');
end

%% Calling the function readDCMfolder
[vol, info] = readDCMfolder();

%% Contrast parameters
global contrastCoronal;
contrastCoronal = [min(vol(:)), max(vol(:))];

%% Display a coronal slice for cropping

hFig = figure('Name', '3D DICOM Viewer with Uniform Size', 'NumberTitle', 'off');
colormap gray;
% Axes for coronal view
hCoronal = axes('Parent', hFig, 'Position', [0.1, 0.2, 0.8, 0.7]);
hCoronalImg = imshow(squeeze(vol(1, :, :))', contrastCoronal, 'Parent', hCoronal);
title(hCoronal, sprintf('Coronal Slice: 1/%d', size(vol, 1)));
axis on;               % Turn on the axis
axis image; 
% Coronal slice slider
uicontrol('Style', 'slider', ...
          'Min', 1, 'Max', size(vol, 1), 'Value', 1, ...
          'SliderStep', [1 / (size(vol, 1) - 1), 10 / (size(vol, 1) - 1)], ...
          'Units', 'normalized', ...
          'Position', [0.4, 0.05, 0.2, 0.05], ...
          'Callback', @(src, ~) updateCoronal(src, vol, hCoronalImg, hCoronal));

%% Function to update coronal slice
function updateCoronal(slider, vol, imgHandle, axHandle)
    % Get the current slice index from the slider
    sliceIndex = round(slider.Value);
    % Update the image data
    set(imgHandle, 'CData', squeeze(vol(sliceIndex, :, :))');
    % Update the title
    title(axHandle, sprintf('Coronal Slice: %d/%d', sliceIndex, size(vol, 1)));
end

figure;
coronal_slice = squeeze(vol(263, :, :))';

imshow(coronal_slice, contrastCoronal);
title('Coronal Slice 263');
axis on;               % Turn on the axis
axis image; 
% Uncomment the following lines for interactive cropping
% cropped_rect = imrect;
% position = wait(cropped_rect);
position=[133.0326   , 0.8739 , 281.3824  , 80.8612];
 % Predefined cropping position
cropped_slice = imcrop(coronal_slice, position);
% resized_slice = imresize(coronal_slice, [size(vol, 1), size(vol, 2)]);
% Get cropping dimensions
cropHeight = round(position(4));
cropWidth = round(position(3));

% Display cropped image
figure;
imshow(cropped_slice, contrastCoronal);
title('Cropped Coronal Slice 263','FontSize',18);
axis on;               % Turn on the axis
axis image; 
%% Segmentation
Igray = mat2gray(cropped_slice);
Th = graythresh(Igray); % Otsu's method4
Ibin = imbinarize(Igray, Th);
Ibin = ~Ibin;

BW = imfill(Ibin, 'holes');
se = strel('disk', 2);
BW_erod = imdilate(BW, se);

%% Step 2: Allow selection of pulmonary branches
figure;
imshow(BW_erod, []);
title('Select Pulmonary Branch Regions (Red)');
axis on;               % Turn on the axis
axis image; 
BW_pulmonary = bwselect;

% Save selected regions as "red regions"
labeledRedRegions = bwlabel(BW_pulmonary);
% labeledRedRegions=imresize(labeledRedRegions, [cropHeight, cropWidth]);
% Overlay the "red regions"
redOverlay = labeloverlay(Igray, labeledRedRegions, 'Colormap', [1 0 0]);
figure;
imshow(redOverlay);
title('Pulmonary Branch Regions Highlighted in Red','FontSize',18);
axis on;               % Turn on the axis
axis image; 

%% Step 3: Selection of lungs
figure;
imshow(BW, []);
title('Select the two lungs and press enter','FontSize',18);
axis on;               % Turn on the axis
axis image; 
BW_selected = bwselect;
BW_selected=imresize(BW_selected, [cropHeight, cropWidth]);
% Calculate region properties
labeledImage0 = bwlabel(BW_selected);
stats = regionprops(labeledImage0, 'Area');
lung_area1 = stats(1).Area;
lung_area2 = stats(2).Area;
figure;
imshow(BW_selected, []);
title('The selected lungs','FontSize',18);
axis on;               % Turn on the axis
axis image; 
%% Create 3D matrix with coronal slices containing the lungs
startSlice = 157;
endSlice = 313;

% Preallocate 3D matrix for cropped coronal slices
croppedCoronalVolume = zeros(cropHeight, cropWidth, endSlice - startSlice +1, 'like', vol);

for sliceIndex = startSlice:endSlice
    coronalSlice = squeeze(vol(sliceIndex, :, :))';
    croppedSlice = imcrop(coronalSlice, position);
    croppedSlice= imresize(croppedSlice, [cropHeight, cropWidth]);
    if size(croppedSlice, 1) ~= size(croppedCoronalVolume, 1) || size(croppedSlice, 2) ~= size(croppedCoronalVolume, 2)
        error('Inconsistent dimensions: Cropped slice is %dx%d but expected %dx%d.', ...
              size(croppedSlice, 1), size(croppedSlice, 2), size(croppedCoronalVolume, 1), size(croppedCoronalVolume, 2));
    end

    croppedCoronalVolume(:, :, sliceIndex - startSlice + 1) = croppedSlice;
end

%% Step 4: Segment lungs in each slice
reference_segmentation = BW_selected;
se = strel('diamond', 1);

refinedVolume = zeros(cropHeight, cropWidth, endSlice - startSlice + 1, 'logical');
Area_lungs = zeros(2, endSlice - startSlice + 1);

for slice_num = startSlice:endSlice
    current_slice = squeeze(vol(slice_num, :, :))';
    cropped_slice = imcrop(current_slice, position);
    cropped_slice=imresize(cropped_slice, [cropHeight, cropWidth]);
    % Binarize the cropped slice
    Igray = mat2gray(cropped_slice);
    Th = graythresh(Igray);
    Ibin = imbinarize(Igray, Th);
    Ibin = ~Ibin;
    I_eroded = imerode(Ibin, se);

    % Refine segmentation
    Iopen = imopen(I_eroded, se);
    labeledImage = bwlabel(Iopen);
    % Resize labeledRedRegions to match cropped slice dimensions
    labeledRedRegions = imresize(labeledRedRegions, [cropHeight, cropWidth], 'nearest');


    refinedMask = false(cropHeight, cropWidth);

    % Process each labeled region
    stats = regionprops(labeledImage, 'Area', 'PixelIdxList');
    for labelIdx = 1:numel(stats)
        region_mask = false(size(labeledImage));
        region_mask(stats(labelIdx).PixelIdxList) = true;

        if stats(labelIdx).Area > 200
            overlap = sum(region_mask(:) & reference_segmentation(:));
            if overlap > 250
                refinedMask = refinedMask | region_mask;
            end
            if (slice_num >= 250) && (slice_num < 279)
                refinedMask = refinedMask & ~labeledRedRegions;
            end
        end
    end

    refinedVolume(:, :, slice_num - startSlice + 1) = refinedMask;
    Area_lungs(1, slice_num - startSlice + 1) = slice_num;
    Area_lungs(2, slice_num - startSlice + 1) = sum(refinedMask(:));
end
%% Visualize all slices with a slider
figure('Name', 'Refined Coronal Slices Viewer', 'NumberTitle', 'off');
colormap gray;

% Display the initial slice
hImg = imshow(refinedVolume(:, :, 1), []);
title(sprintf('Refined Coronal Slice %d', startSlice));
axis on;               % Turn on the axis
axis image; 
% Add a slider for scrolling through the slices
uicontrol('Style', 'slider', ...
          'Min', 1, 'Max', size(refinedVolume, 3), 'Value', 1, ...
          'SliderStep', [1 / (size(refinedVolume, 3) - 1), 10 / (size(refinedVolume, 3) - 1)], ...
          'Units', 'normalized', ...
          'Position', [0.1, 0.02, 0.8, 0.05], ...
          'Callback', @(src, ~) updateSlice(src, hImg, refinedVolume, startSlice));

%% Update function for the slider
function updateSlice(src, hImg, volume, startSlice)
    sliceIndex = round(src.Value);
    set(hImg, 'CData', volume(:, :, sliceIndex));
    title(sprintf('Refined Coronal Slice %d', startSlice + sliceIndex - 1));
end
%% Step5: Area and volume Calculations
%% Calcul des métriques après la boucle
% Dimensions originales
pixel_area = info.PixelSpacing(1) * info.PixelSpacing(2); % mm²/pixel
slice_thickness = info.SliceThickness; % Épaisseur entre les slices (en mm)

% Calcul du volume pulmonaire total (en mm³)
total_volume = sum(Area_lungs(2, :) * pixel_area * slice_thickness);

% Calcul de la hauteur totale des poumons (en cm)
lung_height_mm = size(refinedVolume, 1) * info.PixelSpacing(2);
lung_height_cm = lung_height_mm / 10;

% Calcul de la largeur maximale des poumons (en mm)
lung_widths = zeros(1, size(refinedVolume, 3));
for i = 1:size(refinedVolume, 3)
    [rows, cols] = find(refinedVolume(:, :, i));
    if ~isempty(cols)
        lung_widths(i) = (max(cols) - min(cols)) * info.PixelSpacing(1);
    end
end
max_lung_width_mm = max(lung_widths);


%% Affichage des résultats
% Résumé dans la console
fprintf('Total estimated volume : %.2f mm³ (%.2f l).\n', total_volume, total_volume / 1e6);
fprintf('Lungs height : %.2f cm.\n', lung_height_cm);
fprintf('Max lungs width : %.2f mm.\n', max_lung_width_mm);

% Graphique de la largeur en fonction des slices
figure;
plot(startSlice:endSlice, lung_widths, '-o');
xlabel('Slice number');
ylabel('Width (mm)');
title('Lungs width (Coronal plan)','FontSize',20);
grid on;

% pixel_area = info.PixelSpacing(1) * info.PixelSpacing(2);
% slice_thickness = info.SliceThickness;
% 
% total_volume = sum(Area_lungs(2, :) * pixel_area * slice_thickness);
% lung_height_mm = size(refinedVolume, 1) * info.PixelSpacing(2);
% lung_height_cm = lung_height_mm / 10;


%% Step6 :Volume Viewer
volumeViewer(refinedVolume);

