function [vol, info] = readDCMfolder()
    %% Load DICOM files
    D = dir('*.dcm');
    if isempty(D)
        error('No DICOM files found in the current directory.');
    end

    % Initialize volume and metadata
    im = dicomread(D(1).name);
    vol = zeros(size(im, 1), size(im, 2), numel(D));
    info = dicominfo(D(1).name);

    % Load slices with progress bar
    f = waitbar(0, sprintf('Loading:  %u / %u', 0, numel(D)));
    for ind = 1:numel(D)
        waitbar(ind / numel(D), f, sprintf('Loading:  %u / %u', ind, numel(D)));
        vol(:, :, ind) = dicomread(D(ind).name);
    end
    close(f);


end
