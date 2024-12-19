% Code organization :
%   1. Reading and selecting data
% First, it loads a 3D volume of medical data in DICOM format. A user interface lets you browse the slices and select one for segmentation.
%   2. Segmentation
% The region of interest (ROI) is defined to focus analysis on the lungs and optimize segmentation.
%       2. a) Automatic segmentation :
% Grayscale conversion and use of automatic thresholding (Otsu method) and morphological operations to isolate lungs. The coronal slice area (in mm²) is calculated.
%       2. b) Manual segmentation:
% The user selects the contours of each lung. Masks are created from clicked and interpolated points. The coronal slice area is also calculated.
%   3. Comparison of the two segmentations
% - Dice coefficient: Measures the overall overlap between manual and automatic masks.
% - Average contour distance (ASD): Evaluates average differences between contours.
% - Relative area difference: Compares segmented areas as a percentage.
% - Local curvature analysis: Compares contour curvatures to identify smoothed or diverging areas.

%% Main Function : ctLungSegmentationCoronal
% This function is responsible of code execution.

function ctLungSegmentationCoronal()
    clear all; close all; clc;

    % Reading DICOM files
    [vol, info] = readDCMfolder();

    % Coronal slice selection via interactive interface
    [selected_slice, selected_slice_idx] = visualizeAndSelectCoronalSlice(vol);

    % ROI selection
    [cropped_slice] = cropCoronalSlice(selected_slice);

    % Automatic segmentation
    [auto_mask, coronal_slice_area_auto] = segmentLungAutomatically(cropped_slice, info);

    % Manual segmentation
    [manual_mask, coronal_slice_area_manual] = segmentLungManually(cropped_slice, info);

    % Overlay contours of manual and automatic masks on the cropped image
    overlayContoursCoronalCropped(cropped_slice, manual_mask, auto_mask);

    % Dice coefficient calculation
    dice_score = calculateDice(auto_mask, manual_mask);
    fprintf('Dice coefficient : %.5f\n', dice_score);

    % Calculation of additional comparison variables between segments
    [avgSurfaceDist, areaDiffRel, curvatureDiff] = compareSegmentations(manual_mask, auto_mask, info);

    % Displaying comparison metrics in the console
    fprintf('Average contour distance (ASD) : %.2f mm\n', avgSurfaceDist);
    fprintf('Relative difference in area : %.2f %%\n', areaDiffRel);
    fprintf('Average difference in local curvatures : %.2f\n', curvatureDiff);

    % Results display
    displayResultsCoronal(selected_slice, cropped_slice, auto_mask, manual_mask, ...
        coronal_slice_area_auto, coronal_slice_area_manual, dice_score, ...
        avgSurfaceDist, areaDiffRel, curvatureDiff);
end

%% Function for visualize and select the coronal slice of interest : visualizeAndSelectCoronalSlice
% This function creates an interface for viewing a 3D volume slice by slice, and allows the user to easily select a specific slice for subsequent segmentation.

function [selected_slice, selected_slice_idx] = visualizeAndSelectCoronalSlice(vol)
    contrastCoronal = [min(vol(:)), max(vol(:))];

    num_slices = size(vol, 1);
    current_slice_idx = 1;

    selected_slice = []; % Selected slice
    selected_slice_idx = []; % Index of the Selected slice

    fig = figure('Name', 'Selection of Coronal Slice', 'NumberTitle', 'off', ...
        'Position', [100, 100, 600, 600]); 

    % Interface display
    h_img = imshow(squeeze(vol(current_slice_idx, :, :))', contrastCoronal, 'InitialMagnification', 'fit');
    title(['Slice Coronal: ', num2str(current_slice_idx), '/', num2str(num_slices)]);

    % Create a slider to navigate through slices
    uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', num_slices, 'Value', current_slice_idx, ...
        'SliderStep', [1/(num_slices-1), 10/(num_slices-1)], ...
        'Units', 'normalized', 'Position', [0.2, 0.05, 0.6, 0.03], ...
        'Callback', @sliderCallback);

    % Button to select the slice of interest
    uicontrol('Style', 'pushbutton', 'String', 'Select', ...
        'Units', 'normalized', 'Position', [0.4, 0.9, 0.2, 0.05], ...
        'Callback', @selectSliceCallback);

    setappdata(fig, 'vol', vol);
    setappdata(fig, 'current_slice_idx', current_slice_idx);
    setappdata(fig, 'h_img', h_img);
    setappdata(fig, 'selected', false);

    uiwait(fig); % Wait for user to select slice

    function sliderCallback(src, ~)
        current_slice_idx = round(src.Value);
        setappdata(fig, 'current_slice_idx', current_slice_idx);

        set(h_img, 'CData', squeeze(vol(current_slice_idx, :, :))');
        title(['Slice Coronal: ', num2str(current_slice_idx), '/', num2str(num_slices)]);
    end

    function selectSliceCallback(~, ~)
        selected_slice_idx = getappdata(fig, 'current_slice_idx');
        selected_slice = squeeze(vol(selected_slice_idx, :, :))';

        setappdata(fig, 'selected', true);
        
        close(fig);
    end

    % Closing the figure
    function figureCloseCallback(~, ~)
        if ~getappdata(fig, 'selected')
            selected_slice = []; 
            selected_slice_idx = []; 
        end
        delete(fig);
    end
end


%% Function for selection of the ROI : cropCoronalSlice
% This function crops a specific coronal slice according to predefined dimensions. 

function [cropped_slice] = cropCoronalSlice(slice)
    position = [129.1755, -9.2942, 291.1240, 285.7203]; % Define coordinates of area of interest (rectangle)
    cropped_slice = imcrop(slice, position);

    figure;
    imshow(cropped_slice, []);
    title('Cropped Coronal Slice');
end

%% Function for automatic segmentation : segmentLungAutomatically
% This function performs automatic segmentation. 

function [mask, coronal_slice_area] = segmentLungAutomatically(cropped_slice, info)
    Igray = mat2gray(cropped_slice); % Grayscale conversion

    % Automatic thresholding segmentation
    Th = graythresh(Igray);
    Ibin = imbinarize(Igray, Th);
    Ibin = ~Ibin;
    Ibin = imfill(Ibin, 'holes');

    mask = bwareaopen(Ibin, 500);

    % Calcul of the coronal slice area
    pixel_area = info.PixelSpacing(1) * info.PixelSpacing(2); 
    coronal_slice_area = sum(mask(:)) * pixel_area; % Total lung area in mm²
end

%% Function for manual segmentation : segmentLungManually
% This function performs manual segmentation.

function [mask, coronal_slice_area] = segmentLungManually(cropped_slice, info)
    % Initializing the final mask
    mask = false(size(cropped_slice));
    
    % Selecting points for the first lung
    figure;
    imshow(cropped_slice, []);
    title('Click to select the contours of the first lung (double-click or click “enter” to finish)');
    [x1, y1] = getpts(); 
    poly1 = poly2mask(x1, y1, size(cropped_slice, 1), size(cropped_slice, 2)); 
    
    % Filling holes and cleaning small objects
    poly1 = imfill(poly1, 'holes');   % Fill holes
    poly1 = bwareaopen(poly1, 100);   % Remove small regions with less than 100 pixels

    % Selecting points for the second lung
    figure;
    imshow(cropped_slice, []);
    title('Click to select the contours of the second lung (double-click or click “enter” to finish)');
    [x2, y2] = getpts(); 
    poly2 = poly2mask(x2, y2, size(cropped_slice, 1), size(cropped_slice, 2)); 
    
    % Filling holes and cleaning small objects
    poly2 = imfill(poly2, 'holes');   % Fill holes
    poly2 = bwareaopen(poly2, 100);   % Remove small regions with less than 100 pixels

    % Combine the masks for both lungs
    mask = poly1 | poly2;

    % Calculating the total coronal slice area
    pixel_area = info.PixelSpacing(1) * info.PixelSpacing(2); 
    coronal_slice_area = sum(mask(:)) * pixel_area; % Total lung area in mm²

    % Display final mask superimposed on image
    figure;
    imshow(cropped_slice, []);
    hold on;
    h = imshow(mask);
    set(h, 'AlphaData', 0.5); % Transparency
    title(['Superimposed final mask - Total area: ', num2str(coronal_slice_area), ' mm²']);
end

%% Functions for calculating comparison variables
% This section defines various comparison metrics.

    %% Function to calculate the Dice coefficient : calculateDice
% The Dice coefficient measures the similarity between two sets (or binary masks) by calculating the ratio between the double of their intersection and the sum of their sizes. 
% It varies between 0 (no overlap) and 1 (perfect overlap), and is used to evaluate the quality of the segmentation.

function dice_score = calculateDice(mask1, mask2)
    if size(mask1) ~= size(mask2)
        error('Masks must have the same dimensions to calculate the Dice coefficient');
    end

    % Calcul of the Dice coefficient
    intersection = sum(mask1(:) & mask2(:));
    dice_score = (2 * intersection) / (sum(mask1(:)) + sum(mask2(:)));
end

    %% Function to calculate the comparison metrics : compareSegmentations
% This function calls the other metric calculation functions.

function [avgSurfaceDist, areaDiffRel, curvatureDiff] = compareSegmentations(manualMask, autoMask, info)
    % Extract contours from manual and automatic masks
    manualBoundary = bwboundaries(manualMask);
    autoBoundary = bwboundaries(autoMask);

    % A single contour in each mask
    manualBoundary = manualBoundary{1};
    autoBoundary = autoBoundary{1};

    % Interpolate contours to the same length
    numPoints = 100; 
    manualBoundaryInterp = interp1(1:size(manualBoundary, 1), manualBoundary, linspace(1, size(manualBoundary, 1), numPoints));
    autoBoundaryInterp = interp1(1:size(autoBoundary, 1), autoBoundary, linspace(1, size(autoBoundary, 1), numPoints));

    % Calculation of average contour distance (ASD)
    avgSurfaceDist = computeAverageSurfaceDist(manualBoundaryInterp, autoBoundaryInterp);

    % Comparison of segmented areas
    pixelArea = info.PixelSpacing(1) * info.PixelSpacing(2); 
    manualArea = sum(manualMask(:)) * pixelArea;
    autoArea = sum(autoMask(:)) * pixelArea;
    areaDiffRel = abs(manualArea - autoArea) / manualArea * 100; % Relative difference in %

    % Local curvature calculation
    manualCurvature = computeLocalCurvature(manualBoundaryInterp);
    autoCurvature = computeLocalCurvature(autoBoundaryInterp);

    % Curvature comparison (average difference in local curvatures)
    curvatureDiff = mean(abs(manualCurvature - autoCurvature));

    % Results display
    figure;
    subplot(1, 2, 1);
    plot(manualCurvature, 'LineWidth', 1.5);
    title('Local curvature (manual)');
    xlabel('Point on contour');
    ylabel('Curvature');
    grid on;

    subplot(1, 2, 2);
    plot(autoCurvature, 'LineWidth', 1.5);
    title('Local curvature (automatic)');
    xlabel('Point on contour');
    ylabel('Curvature');
    grid on;
end

    %% Function to calculate the average surface distance (ASD) : computeAverageSurfaceDist
% The Average Surface Distance (ASD) measures the average distance between the surfaces of the two segmented regions. 
function avgSurfaceDist = computeAverageSurfaceDist(boundary1, boundary2) % Calculates minimum distances between points on two contours 
    d1 = pdist2(boundary1, boundary2, 'euclidean'); 
    d2 = pdist2(boundary2, boundary1, 'euclidean');

    % Average minimum distances
    avgDist1 = mean(min(d1, [], 2));
    avgDist2 = mean(min(d2, [], 2));
    
    % Overall average (symmetrical)
    avgSurfaceDist = (avgDist1 + avgDist2) / 2;
end

    %% Function to compare curvature : computeLocalCurvature
% Curvature comparison quantifies how a contour curves at each point, using first derivatives (slopes) and second derivatives of coordinates to measure local variations in contour direction and shape.

function curvature = computeLocalCurvature(boundary)
    x = boundary(:, 2); 
    y = boundary(:, 1);

    % Calculation of first derivatives (centered approximations)
    dx = gradient(x);
    dy = gradient(y);

    % Calculating second derivatives
    d2x = gradient(dx);
    d2y = gradient(dy);

    % Curvature formula : κ = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
    numerator = abs(dx .* d2y - dy .* d2x);
    denominator = (dx.^2 + dy.^2).^(3/2);

    curvature = numerator ./ (denominator + eps); % Avoid division by zero

    curvature(isnan(curvature)) = 0; % Fill NaN with 0 if isolated points are present
end

%% Function to display results : displayResultsCoronal
% This function is used to visually display and summarize the results of coronal segmentation analysis.

function displayResultsCoronal(slice, cropped_slice, auto_mask, manual_mask, ...
    auto_area, manual_area, dice, avgSurfaceDist, areaDiffRel, curvatureDiff)

    figure;

    % Display of the original slice
    subplot(2, 3, 1);
    imshow(slice, []);
    title('Original Coronal Slice');

    % Displaying the cropped slice
    subplot(2, 3, 2);
    imshow(cropped_slice, []);
    title('Cropped Slice');

    % Automatic segmentation
    subplot(2, 3, 4);
    imshow(auto_mask, []);
    title(['Auto : ', num2str(auto_area), ' mm²']);

    % Manual segmentation
    subplot(2, 3, 5);
    imshow(manual_mask, []);
    title(['Manual : ', num2str(manual_area), ' mm²']);

    % Displaying metrics in a subplot
    subplot(2, 3, 6);
    metricsText = {
        ['Dice : ', num2str(dice)], ...
        ['ASD : ', num2str(avgSurfaceDist), ' mm'], ...
        ['Relative error : ', num2str(areaDiffRel), ' %'], ...
        ['Curvature : ', num2str(curvatureDiff), ' (average)']};
    text(0.5, 0.5, metricsText, 'FontSize', 12, 'HorizontalAlignment', 'center');
    axis off;
end

%% Function to overlay contours on the cropped coronal image
function overlayContoursCoronalCropped(cropped_slice, manual_mask, auto_mask)
    % Extract boundaries for manual and automatic masks
    manualBoundary = bwboundaries(manual_mask);
    autoBoundary = bwboundaries(auto_mask);
    
    % Display the cropped slice as background
    figure;
    imshow(cropped_slice, []); % Use the cropped image as background
    hold on;
    
    % Plot manual mask contours in blue
    for k = 1:length(manualBoundary)
        boundary = manualBoundary{k};
        plot(boundary(:, 2), boundary(:, 1), 'b-', 'LineWidth', 1.5); % Solid blue line
    end
    
    % Plot automatic mask contours in red
    for k = 1:length(autoBoundary)
        boundary = autoBoundary{k};
        plot(boundary(:, 2), boundary(:, 1), 'r--', 'LineWidth', 1.5); % Red dashed line
    end
    
    % Fix the legend with correct colors and styles
    h1 = plot(NaN, NaN, 'b-', 'LineWidth', 1.5); % Dummy for manual mask
    h2 = plot(NaN, NaN, 'r--', 'LineWidth', 1.5); % Dummy for automatic mask
    
    legend([h1, h2], {'Manual Mask', 'Automatic Mask'}, 'Location', 'southoutside', ...
           'Orientation', 'horizontal', 'FontSize', 10);
    title('Overlay of Manual (Blue) and Automatic (Red) Contours - Cropped Coronal');
    hold off;
end







