% This only thng to change for using the code is the path to the files
% (2 last lines of the code and the first function)

%% Function ctLungSegmentation
%This function performs lung segmentaton on a CT scan by reading DICOM files, allowing the user to select a specific slice, and returning the selected slice and its index.
function [selected_slice, selected_slice_idx] = ctLungSegmentation()
    clear all; close all; clc;
    
    %Folder path to be modified
    folder_path = '/Users/maelebedel/Desktop/Politecnico/Cours/Biomedical_signal_and_Images_processing/Biomedical_Images_Processing/Projet/Projet_BIP/CT';
    
    %Read DICOM files by calling the readDCMfolder function
    [image_ct, info] = readDCMfolder();
    num_slices = size(image_ct, 3); 
    %Call the visualizeAndSelectSlice function, which allows the user to select one
    [selected_slice, selected_slice_idx] = visualizeAndSelectSlice(image_ct);
end

%% Function visualizeAndSelectSlice
%This function creates an interactive viewer that allows the user to navigate through CT slices using a slider and select a specific slice for further processing.
function [selected_slice, selected_slice_idx] = visualizeAndSelectSlice(image_ct)
    
    %Create figure with slider to navigate through slices
    num_slices = size(image_ct, 3);
    current_slice_idx = 1;
    
    %Create figure
    fig = figure('Name', 'CT Slice Viewer', 'NumberTitle', 'off', ...
        'CloseRequestFcn', @figureCloseCallback);
    
    %Display initial slice
    h_img = imshow(image_ct(:,:,current_slice_idx), [], 'InitialMagnification', 'fit');
    title(['Slice ', num2str(current_slice_idx), '/', num2str(num_slices)]);
    
    %Add slider
    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', num_slices, ...
        'Value', current_slice_idx, ...
        'Units', 'normalized', ...
        'Position', [0.2 0.02 0.6 0.03], ...
        'Callback', @sliderCallback);
    
    %Add select button
    select_btn = uicontrol('Style', 'pushbutton', ...
        'String', 'Select This Slice', ...
        'Units', 'normalized', ...
        'Position', [0.4 0.95 0.2 0.05], ...
        'Callback', @selectSliceCallback);
    
    %Store data in figure's application data
    setappdata(fig, 'image_ct', image_ct);
    setappdata(fig, 'current_slice_idx', current_slice_idx);
    setappdata(fig, 'h_img', h_img);
    setappdata(fig, 'selected', false);
    
    %Slider callback function
    function sliderCallback(src, ~)
        %Get current figure
        fig = ancestor(src, 'figure');
        
        %Update slice index
        current_slice_idx = round(src.Value);
        setappdata(fig, 'current_slice_idx', current_slice_idx);
        
        %Get image_ct and h_img
        image_ct = getappdata(fig, 'image_ct');
        h_img = getappdata(fig, 'h_img');
        
        %Update image
        set(h_img, 'CData', image_ct(:,:,current_slice_idx));
        title(['Slice ', num2str(current_slice_idx), '/', num2str(size(image_ct, 3))]);
    end
    
    %Select slice callback
    function selectSliceCallback(src, ~)
        %Get current figure
        fig = ancestor(src, 'figure');
        
        %Get current slice index and image_ct
        current_slice_idx = getappdata(fig, 'current_slice_idx');
        image_ct = getappdata(fig, 'image_ct');
        
        %Mark as selected and close
        setappdata(fig, 'selected', true);
        close(fig);
    end
    
    %Figure close callback
    function figureCloseCallback(src, ~)
        %Check if a slice was selected
        if ~getappdata(src, 'selected')
            %If not selected, cancel the operation
            selected_slice = [];
            selected_slice_idx = [];
            delete(src);
            return;
        end
        
        %Get selected slice information
        current_slice_idx = getappdata(src, 'current_slice_idx');
        image_ct = getappdata(src, 'image_ct');
        
        %Set output variables
        selected_slice = image_ct(:,:,current_slice_idx);
        selected_slice_idx = current_slice_idx;
        
        %Delete the figure
        delete(src);
    end
    
    % ait for user to select a slice
    uiwait(fig);
end


%% Function to add some kind of noise
%This function adds noise to an input image, supporting 'pepper', 'gaussian', or no noise types, with parameters specifying the noise characteristics.
function noisy_image = addNoise(input_image, noise_type, noise_params)
    
    img = double(input_image);

    switch noise_type
        case 'pepper'
            %Add pepper noise
            noise_density = noise_params;  %Density (0 to 1)
            rand_matrix = rand(size(img));
            noisy_image = img;
            noisy_image(rand_matrix < noise_density) = 0;

        case 'gaussian'
            %Add Gaussian noise
            mean = noise_params(1);    %Mean 
            std_dev = noise_params(2); 
            noisy_image = img + std_dev * randn(size(img)) + mean;

        case 'none'
            %No noise, return the original image
            noisy_image = img;

        otherwise
            error('Unknown noise type. Choose ''pepper'', ''gaussian'', or ''none''.');
    end
end

%% Function selection_roi_image1
%This function selects a predefined rectangular region of interest (ROI) from an input image, adds specified noise, and returns the cropped and normalized ROI.

function [rect, cropped_slice] = selection_roi_image1(slice, noise_type, parameter)
    slice_double1 = double(slice);
    slice_double = addNoise(slice_double1, noise_type, parameter);
    rect =[124.0000;  127.0000;  300.0000;  216.0000];
    
    %Extract the region of interest (ROI) from the image
    cropped_slice1 = slice_double(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1);
    cropped_slice = mat2gray(cropped_slice1);  
end

%% Function region_interest
%This function allows the user to select two regions of interest in a binary image (e.g., lungs) and creates a mask that combines these regions.
function [mask] = region_interet(slice)

    %Apply a threshold to binarize the image
    Th = graythresh(slice);
    bSlice = imbinarize(slice, Th);
    bSlice = ~bSlice;

    %Label the segmented regions
    labeled = bwlabel(bSlice, 8);   
    colored_labels = label2rgb(labeled, 'hsv', 'k', 'shuffle');

    %Display the image for the user to select two points
    figure;
    imshow(colored_labels, []);
    title('Click on two regions representing the lungs');
    
    %Ask the user to select two points (one for each lung)
    [x, y] = ginput(2);  
    if numel(x) ~= 2 || numel(y) ~= 2
        disp('You must select two points.');
        mask = [];
        return;
    end
    
    %Create a mask combining the regions containing the two selected points
    mask = false(size(bSlice));
    for i = 1:2
        mask = mask | bwselect(labeled, x(i), y(i), 8);
    end
    
    figure;
    imshow(mask, []);
    title('Mask defined by the two selected regions');
end

%% Function applySegmentationToAllSlices
%This function processes all CT slices from a specified folder, applies noise, segments lung regions, calculates cross-sectional areas, and estimates the total lung volume and height based on the segmented slices.



function applySegmentationToAllSlices(folder_path)
    
    %Call the lung CT segmentation function to choose the first one
    [selected_slice, selected_slice_idx] = ctLungSegmentation();
    
    %Ask the user to select noise type and parameters
    noise_types = {'pepper', 'gaussian', 'none'};
    [noise_type_idx, ok] = listdlg('PromptString', 'Select the noise type:', ...
        'SelectionMode', 'single', 'ListString', noise_types);
    
    if ~ok
        disp('Operation canceled by user.');
        return;
    end
    
    noise_type = noise_types{noise_type_idx};
    
    %Get noise parameters from user
    if strcmp(noise_type, 'pepper')
        noise_param = inputdlg('Enter the noise density (0 to 1):', 'Pepper Noise Density', [1 50], {'0.2'});
        noise_params = str2double(noise_param{1});
    elseif strcmp(noise_type, 'gaussian')
        noise_params = inputdlg({'Enter the mean of Gaussian noise:', 'Enter the standard deviation of Gaussian noise:'}, ...
            'Gaussian Noise Parameters', [1 50], {'0', '500'});
        noise_params = str2double(noise_params);
    else
        noise_params = []; % No noise parameters for 'none'
    end

    %Call te ROI fonction to define the rectangle crop
    [rect, cropped_slice] = selection_roi_image1(selected_slice,noise_type, noise_params);
    
    %Creation of a mask based on the ROI
    initial_mask = region_interet(cropped_slice);

    %Read the DICOM files from the folder
    [image_ct, info] = readDCMfolder();  
    num_slices = min(size(image_ct, 3), 63);  % Limit the number of slices to 66 because we don't have lungs after
    y_limit = 100;
    
    %Array to store cross-sectional areas
    slice_areas = [];
    
    %Initialize a 3D volume for visualization
    slice_spacing = 4; % Espace entre les tranches en mm
    num_interpolated_slices = round(slice_spacing / info.SliceThickness) - 1;
    volume = false(size(image_ct, 1), size(image_ct, 2), num_slices + num_interpolated_slices * (num_slices - 1));

    %Loop through all slices
    for idx = 1:num_slices

        %Get the current slice
        current_slice = image_ct(:,:,idx);
        current_slice_double1 = double(current_slice);
        current_slice_double = addNoise(current_slice_double1, noise_type, noise_params);
        
        %Extract the region of interest (ROI) from the image based oh the rectangle of the first image
        current_cropped_slice = current_slice_double(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1);
        current_cropped_slice = mat2gray(current_cropped_slice);
        
        %Apply a threshold to binarize the image
        Th = graythresh(current_cropped_slice);
        bSlice = imbinarize(current_cropped_slice, Th);
        bSlice = ~bSlice;
        bSlice = imfill(bSlice, 'holes');  % Fill holes
        bSlice = bwareaopen(bSlice, 70);  % Remove small noisy regions

        if 2< idx < 54
            %Label the segmented regions
            labeled = bwlabel(bSlice, 8);

            %Retain regions that intersect with the initial mask
            selected_regions = false(size(bSlice));
            unique_labels = unique(labeled(initial_mask));
            for label = unique_labels(:)'
                if label > 0 % Avoid label 0 (background)
                    selected_regions = selected_regions | (labeled == label);
                end
            end

            %Identify regions in selected_regions
            labeled_selected = bwlabel(selected_regions, 8);
            stats = regionprops(labeled_selected, 'Area');  
            areas = [stats.Area];  

            %index of the 2 bigger area (lungs)
            [~, sorted_indices] = sort(areas, 'descend');
            if numel(sorted_indices) > 2
                largest_indices = sorted_indices(1:2); 
            else
                largest_indices = sorted_indices; 
            end
            selected_regions = false(size(selected_regions));
            for i = 1:numel(largest_indices)
                selected_regions = selected_regions | (labeled_selected == largest_indices(i));
            end
            
        elseif idx >= 54           

            %Apply the mask_no_top to exclude the top part
            bSlice(1:y_limit, :) = false;  % Exclude the top part

            %Label the segmented regions
            labeled = bwlabel(bSlice, 8);

            %Retain regions that intersect with the initial mask and the no_top_mask
            selected_regions = false(size(bSlice));
            unique_labels = unique(labeled(initial_mask)); 
            for label = unique_labels(:)'
                if label > 0 
                    selected_regions = selected_regions | (labeled == label);
                end
            end
        end

        %Resize selected_regions to match the original image size for the calculation of the volume
        selected_regions_resized = false(size(current_slice));
        [rows, cols] = size(selected_regions);
        selected_regions_resized(1:rows, 1:cols) = selected_regions;

        %Add the selected regions to the 3D volume
        volume(:,:,1 + (idx-1) * (1 + num_interpolated_slices)) = selected_regions_resized;
        if idx < num_slices
            for k = 1:num_interpolated_slices
                volume(:,:,1 + (idx-1) * (1 + num_interpolated_slices) + k) = false(size(selected_regions_resized));
            end
        end

        %Calculate cross-sectional area
        cross_sectional_area = sum(selected_regions(:));
        slice_areas = [slice_areas; idx, cross_sectional_area]; % Append results
    end

    %Find the slice with the largest cross-sectional area
    [max_area, max_area_idx] = max(slice_areas(:, 2)); %Maximum area and its index
    slice_with_max_area = slice_areas(max_area_idx, 1); % Slice number corresponding to the maximum area
    

    %Calculate total lung volume
    slice_thickness = info.SliceThickness; %Thickness between slices (in mm)
    pixel_area = info.PixelSpacing(1) * info.PixelSpacing(2); %Area of a pixel (in mm²): converts each pixel to its equivalent in mm² (calculated from info.PixelSpacing)
    total_volume = sum(slice_areas(:, 2) * pixel_area * slice_thickness); % Total volume (in mm³): slice_area(:,2) represents the number of segmented pixels converted to mm² using pixel_area, and slice_thickness is the thickness between slices in mm

    
    %Calculate total lung height
    lung_height_mm = num_slices * slice_thickness; % Height in mm
    lung_height_cm = lung_height_mm / 10; % Conversion to cm

    %Display results in console
    fprintf('Slices containing the lungs and their cross-sectional areas :\n');
    fprintf('Slice Index    Cross-sectional Area (mm²)\n');
    fprintf('------------------------------------------\n');
    for i = 1:size(slice_areas, 1)
        fprintf('%8d %25d\n', slice_areas(i, 1), slice_areas(i, 2));
    end
    fprintf('\nThe slice with the largest cross-sectional area is slice %d with an area of %.2f mm².\n', slice_with_max_area, max_area * pixel_area);
    fprintf('\nEstimated total lung height : %.2f cm.\n', lung_height_cm);
    fprintf('\nEstimated total lung volume : %.2f mm³ (%.2f liters).\n', total_volume, total_volume / 1e6);

    %Plot the cross-sectional areas
    figure;
    plot(slice_areas(:, 1), slice_areas(:, 2), '-o');
    xlabel('Slice Index');
    ylabel('Cross-sectional Area (mm²)');
    title('Cross-sectional Areas of Lung Slices');
    grid on;

    %3D Visualization of the segmented volume
    volumeViewer(volume); 
end

folder_path = 'C:\Users\Farah\Desktop\Cours_Polimi\Signal_Processing\Projet2\CT';
applySegmentationToAllSlices(folder_path);
