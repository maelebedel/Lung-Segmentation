% Code organization :
%   1. Reading and selecting data
% First, it loads a 3D volume of medical data in DICOM format. A user interface lets you browse the slices and select one for segmentation.
%   2. Segmentation
% The region of interest (ROI) is defined to focus analysis on the lungs and optimize segmentation.
%       2. a) Automatic segmentation :
% Grayscale conversion and use of automatic thresholding (Otsu method) and morphological operations to isolate lungs. The cross sectionnal area (in mm²) is calculated.
%       2. b) Manual segmentation:
% The user selects the contours of each lung. Masks are created from clicked and interpolated points. The cross sectionnal area is also calculated.
%   3. Comparison of the two segmentations
% - Dice coefficient: Measures the overall overlap between manual and automatic masks.
% - Average contour distance (ASD): Evaluates average differences between contours.
% - Relative area difference: Compares segmented areas as a percentage.
% - Local curvature analysis: Compares contour curvatures to identify smoothed or diverging areas.

%% Main Function : ctLungSegmentationAxial
function ctLungSegmentationAxial()
    clear all; close all; clc;

    % Reading DICOM files
    [image_ct, info] = readDCMfolder();

    % Axial slice selection via interactive interface
    [selected_slice, selected_slice_idx] = visualizeAndSelectSlice(image_ct);

    % ROI selection
    [cropped_slice] = selectROIFixed(selected_slice);

    % Automatic segmentation
    [auto_mask, cross_sectional_area_auto] = segmentLungAutomatically(cropped_slice, info);

    % Manual segmentation
    [manual_mask, cross_sectional_area_manual] = segmentLungManually(cropped_slice, info);
    
    % Overlay contours of manual and automatic masks
    overlayContours(cropped_slice, manual_mask, auto_mask);

    % Comparison of manual and automatic segmentations
    [avgSurfaceDist, areaDiffRel, curvatureDiff, manualCurvature, autoCurvature] = compareSegmentations(manual_mask, auto_mask, info);

    % Displaying comparison metrics in the console
    disp(['Cross-sectional area (automatic) : ', num2str(cross_sectional_area_auto), ' mm²']);
    disp(['Cross-sectional area (manual) : ', num2str(cross_sectional_area_manual), ' mm²']);
    fprintf('Dice coefficient : %.5f\n', calculateDice(manual_mask, auto_mask));
    fprintf('Average contour distance (ASD) : %.2f mm\n', avgSurfaceDist);
    fprintf('Relative difference in area : %.2f %%\n', areaDiffRel);
    fprintf('Average difference in local curvatures : %.2f\n', curvatureDiff);

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

    displayResultsAxial(selected_slice, cropped_slice, auto_mask, manual_mask, ...
        cross_sectional_area_auto, cross_sectional_area_manual, calculateDice(manual_mask, auto_mask), ...
        avgSurfaceDist, areaDiffRel, curvatureDiff);
end

%% Manual VS automatic comparison function : compareSegmentations
% This function calls the other metric calculation functions.

function [avgSurfaceDist, areaDiffRel, curvatureDiff, manualCurvature, autoCurvature] = compareSegmentations(manualMask, autoMask, info)
    % Extract mask contours
    manualBoundary = bwboundaries(manualMask);
    autoBoundary = bwboundaries(autoMask);

    % Interpolate contours to have a fixed number of points
    numPoints = 100;
    manualBoundaryInterp = interp1(1:size(manualBoundary{1}, 1), manualBoundary{1}, linspace(1, size(manualBoundary{1}, 1), numPoints));
    autoBoundaryInterp = interp1(1:size(autoBoundary{1}, 1), autoBoundary{1}, linspace(1, size(autoBoundary{1}, 1), numPoints));

    % ASD Calculation
    avgSurfaceDist = computeAverageSurfaceDist(manualBoundaryInterp, autoBoundaryInterp);

    % Calculation of area differences
    pixelArea = info.PixelSpacing(1) * info.PixelSpacing(2); 
    manualArea = sum(manualMask(:)) * pixelArea;
    autoArea = sum(autoMask(:)) * pixelArea;
    areaDiffRel = abs(manualArea - autoArea) / manualArea * 100;

    % Calculation of local curvatures
    manualCurvature = computeLocalCurvature(manualBoundaryInterp);
    autoCurvature = computeLocalCurvature(autoBoundaryInterp);
    curvatureDiff = mean(abs(manualCurvature - autoCurvature));
end

%% Function to display results : displayResultsAxial
% This function is used to visually display and summarize the results of axial segmentation analysis.

function displayResultsAxial(slice, cropped_slice, auto_mask, manual_mask, ...
    auto_area, manual_area, dice, avgSurfaceDist, areaDiffRel, curvatureDiff)

    figure;

    % Display of the original slice
    subplot(2, 3, 1);
    imshow(slice, []);
    title('Original Axial Slice');

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

%% Function to calculate the average surface distance (ASD) : computeAverageSurfaceDist
% The Average Surface Distance (ASD) measures the average distance between the surfaces of the two segmented regions. 

function avgSurfaceDist = computeAverageSurfaceDist(boundary1, boundary2)
    d1 = pdist2(boundary1, boundary2, 'euclidean');
    d2 = pdist2(boundary2, boundary1, 'euclidean');
    avgSurfaceDist = (mean(min(d1, [], 2)) + mean(min(d2, [], 2))) / 2;
end

%% Function to compare curvature : computeLocalCurvature
% Curvature comparison quantifies how a contour curves at each point, using first derivatives (slopes) and second derivatives of coordinates to measure local variations in contour direction and shape.

function curvature = computeLocalCurvature(boundary)
    x = boundary(:, 2);
    y = boundary(:, 1);
    dx = gradient(x);
    dy = gradient(y);
    d2x = gradient(dx);
    d2y = gradient(dy);
    numerator = abs(dx .* d2y - dy .* d2x);
    denominator = (dx.^2 + dy.^2).^(3/2);
    curvature = numerator ./ (denominator + eps);
    curvature(isnan(curvature)) = 0;
end

%% Function for visualize and select the axial slice of interest : visualizeAndSelectSlice
% This function creates an interface for viewing a 3D volume slice by slice, and allows the user to easily select a specific slice for subsequent segmentation.

function [selected_slice, selected_slice_idx] = visualizeAndSelectSlice(image_ct)
    num_slices = size(image_ct, 3);
    current_slice_idx = 1;

    % Interface display
    fig = figure('Name', 'Sélection de Slice', 'NumberTitle', 'off', ...
        'CloseRequestFcn', @figureCloseCallback);
    h_img = imshow(image_ct(:, :, current_slice_idx), [], 'InitialMagnification', 'fit');
    title(['Slice ', num2str(current_slice_idx), '/', num2str(num_slices)]);
    
    uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', num_slices, 'Value', current_slice_idx, ...
        'Units', 'normalized', 'Position', [0.2 0.02 0.6 0.03], ...
        'Callback', @sliderCallback);
    
    uicontrol('Style', 'pushbutton', 'String', 'Sélectionner', ...
        'Units', 'normalized', 'Position', [0.4 0.95 0.2 0.05], ...
        'Callback', @selectSliceCallback);
    
    setappdata(fig, 'image_ct', image_ct);
    setappdata(fig, 'current_slice_idx', current_slice_idx);
    setappdata(fig, 'h_img', h_img);
    setappdata(fig, 'selected', false);
    
    uiwait(fig);
    
    function sliderCallback(src, ~)
        current_slice_idx = round(src.Value);
        setappdata(fig, 'current_slice_idx', current_slice_idx);
        set(h_img, 'CData', image_ct(:, :, current_slice_idx));
        title(['Slice ', num2str(current_slice_idx), '/', num2str(num_slices)]);
    end

    function selectSliceCallback(~, ~)
        setappdata(fig, 'selected', true);
        close(fig);
    end

    function figureCloseCallback(~, ~)
        if getappdata(fig, 'selected')
            selected_slice_idx = getappdata(fig, 'current_slice_idx');
            selected_slice = image_ct(:, :, selected_slice_idx);
        else
            selected_slice_idx = [];
            selected_slice = [];
        end
        delete(fig);
    end
end

%% Function for selection of the ROI : selectROIFixed
% This function crops a specific axial slice according to predefined dimensions. 

function [cropped_slice] = selectROIFixed(slice)
    rect = [124, 127, 300, 216]; % [x, y, width, height]

    slice_double = double(slice);
    cropped_slice1 = slice_double(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1);

    cropped_slice = mat2gray(cropped_slice1);
end

%% Function for automatic segmentation : segmentLungAutomatically
% This function performs automatic segmentation. 

function [mask, cross_sectional_area] = segmentLungAutomatically(cropped_slice, info)
    slice_double = double(cropped_slice);

    % Automatic segmentation with Otsu thresholding
    Th = graythresh(slice_double);
    bSlice = imbinarize(slice_double, Th);
    bSlice = ~bSlice;  % Inversion so that the lungs are white
    bSlice = imfill(bSlice, 'holes');  % Fill the holes
    bSlice = bwareaopen(bSlice, 70);  % Removes small spurious regions

    % Labeling of connected regions
    labeled = bwlabel(bSlice, 8);
    colored_labels = label2rgb(labeled, 'hsv', 'k', 'shuffle');  

    % Interactive display to select lungs
    figure;
    imshow(colored_labels, []);
    title('Click on two regions representing the lungs');
    [x, y] = ginput(2); 

    % Creating the combined mask for the two selected regions
    mask = false(size(bSlice));
    for i = 1:2
        mask = mask | bwselect(labeled, round(x(i)), round(y(i)), 8);
    end

    % Calculation of cross-sectional area
    pixel_area = info.PixelSpacing(1) * info.PixelSpacing(2); % Area of ​​a pixel in mm²
    cross_sectional_area = sum(mask(:)) * pixel_area; % Total lung area in mm²
end

%% Function for manual segmentation : segmentLungManually
% This function performs manual segmentation.

function [mask, cross_sectional_area] = segmentLungManually(cropped_slice, info)
    % Initialize the final mask
    mask = false(size(cropped_slice));
    
    % Selecting points for the first lung
    figure;
    imshow(cropped_slice, []);
    title('Click to select the outlines of the first lung (double-click or click enter to finish)');
    [x1, y1] = getpts(); % Retrieve the clicked points for the first lung
    % Ensure the contour is closed
    x1 = [x1; x1(1)];
    y1 = [y1; y1(1)];
    poly1 = poly2mask(x1, y1, size(cropped_slice, 1), size(cropped_slice, 2));

    % Selecting points for the second lung
    figure;
    imshow(cropped_slice, []);
    title('Click to select the outlines of the second lung (double-click or click enter to finish)');
    [x2, y2] = getpts(); % Retrieve the clicked points for the second lung
    % Ensure the contour is closed
    x2 = [x2; x2(1)];
    y2 = [y2; y2(1)];
    poly2 = poly2mask(x2, y2, size(cropped_slice, 1), size(cropped_slice, 2));

    % Combine the two masks
    mask = poly1 | poly2;

    % Calculating the total cross-sectional area
    pixel_area = info.PixelSpacing(1) * info.PixelSpacing(2);  % Area of a pixel in mm²
    cross_sectional_area = sum(mask(:)) * pixel_area; % Total lung area in mm²

    % Display the final mask superimposed on the cropped slice
    figure;
    imshow(cropped_slice, []);
    hold on;
    h = imshow(mask);
    set(h, 'AlphaData', 0.5); 
    title(['Superimposed final mask - Total area: ', num2str(cross_sectional_area), ' mm²']);
end

%% Function to calculate the Dice coefficient : calculateDice
% The Dice coefficient measures the similarity between two sets (or binary masks) by calculating the ratio between the double of their intersection and the sum of their sizes. 
% It varies between 0 (no overlap) and 1 (perfect overlap), and is used to evaluate the quality of the segmentation.

function dice_score = calculateDice(mask1, mask2)
    % Check that the two masks are the same size
    if size(mask1) ~= size(mask2)
        error('Masks must have the same dimensions to calculate the Dice coefficient');
    end

    % Calcul of the Dice coefficient
    intersection = sum(mask1(:) & mask2(:)); % Pixels common to both masks
    dice_score = (2 * intersection) / (sum(mask1(:)) + sum(mask2(:))); % Dice coefficient formula
end

%% Function to overlay contours of manual and automatic masks
function overlayContours(cropped_slice, manual_mask, auto_mask)
    % Extract boundaries for manual and automatic masks
    manualBoundary = bwboundaries(manual_mask);
    autoBoundary = bwboundaries(auto_mask);
    
    % Display the original cropped slice as background
    figure;
    imshow(cropped_slice, []);
    hold on;
    
    % Plot manual mask contours in blue
    for k = 1:length(manualBoundary)
        boundary = manualBoundary{k};
        plot(boundary(:, 2), boundary(:, 1), 'b-', 'LineWidth', 1.5); % Blue solid line
    end
    
    % Plot automatic mask contours in red
    for k = 1:length(autoBoundary)
        boundary = autoBoundary{k};
        plot(boundary(:, 2), boundary(:, 1), 'r--', 'LineWidth', 1.5); % Red dashed line
    end
    
    % Create explicit handles for the legend
    h1 = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);  % Blue solid line for manual mask
    h2 = plot(NaN, NaN, 'r--', 'LineWidth', 1.5); % Red dashed line for automatic mask

    % Add legend with correct handles
    legend([h1, h2], {'Manual Mask', 'Automatic Mask'}, 'Location', 'northeast');

    title('Overlay of Manual (Blue) and Automatic (Red) Contours');
    hold off;
end




