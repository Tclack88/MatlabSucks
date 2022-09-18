function [t, x, px, dt]=ColourTracking(vidObj)
    
    DEBUG_MASKS = false;
    
    nFrames=get(vidObj, 'NumberOfFrames');
    framerate = get(vidObj, 'FrameRate');
    
    dt = 1/framerate;
    t = 0:dt:(dt*(nFrames-1));
    
    % Read a single frame from the middle of the video and ask the user to pick
    % the colors for the controller and the plans
    [hueControl, huePlant] = promptForMaskCalibrationValues(read(vidObj, floor(nFrames/2)));
    fprintf('Hue for controller: %0.8f, Hue for Plant: %0.8f\n', hueControl, huePlant);

    
    if ((exist('DEBUG_MASKS', 'var') == 1) && DEBUG_MASKS)
        hFigCurrentFrame = figure();
        axCurrentRGB = subplot(1, 3, 1, 'Parent', hFigCurrentFrame);
        hImgCurrentRGB = imshow([], 'Parent', axCurrentRGB);
        title(axCurrentRGB, 'Current RGB Frame');
        
        axCurrentControlMask = subplot(1, 3, 2, 'Parent', hFigCurrentFrame); hold(axCurrentControlMask, 'on');
        hImgCurrentControlMask = imshow([], 'Parent', axCurrentControlMask);
        hImgCurrentControlCenter = plot(axCurrentControlMask, nan, nan, 'gp');
        title(axCurrentControlMask, 'Controller Mask');
        
        axCurrentPlantMask = subplot(1, 3, 3, 'Parent', hFigCurrentFrame); hold(axCurrentPlantMask, 'on');
        hImgCurrentPlantMask = imshow([], 'Parent', axCurrentPlantMask);
        hImgCurrentPlantCenter = plot(axCurrentPlantMask, nan, nan, 'gp');
        title(axCurrentPlantMask, 'Plant Mask');
    end
    
    x = zeros(1, nFrames);
    px = zeros(1, nFrames);
    
    lastBlobControl = [];
    lastBlobPlant = [];
    
    for iFrame=1:nFrames
        fprintf('Processing Frame %u of %u (%4.1f%%)\n', iFrame, nFrames, (iFrame*100)/nFrames);
        
        imgRGB = read(vidObj, iFrame); % Acquire single frame
        imgHSV = rgb2hsv(imgRGB);
        
        %set(hImgCurrentRGB, 'CData', imgRGB);
        
        maskControl = getHueMaskHSV(imgHSV, hueControl);
        maskControl = FillMaskHoles(maskControl);
        blobControl = GetBiggestMaskBlob(maskControl);
        if (isempty(blobControl))
            assert(~isempty(lastBlobControl), 'Cannot accept ''no blob detected'' on first frame');
            blobControl = lastBlobControl;
        else
            lastBlobControl = blobControl;
        end
        
        maskPlant = getHueMaskHSV(imgHSV, huePlant);
        maskPlant = FillMaskHoles(maskPlant);
        blobPlant = GetBiggestMaskBlob(maskPlant);
        if (isempty(blobPlant))
            assert(~isempty(lastBlobPlant), 'Cannot accept ''no blob detected'' on first frame');
            blobPlant = lastBlobPlant;
        else
            lastBlobPlant = blobPlant;
        end
        
        
        x(iFrame) = blobControl(2);
        px(iFrame) = blobPlant(2);
        
        if ((exist('DEBUG_MASKS', 'var') == 1) && DEBUG_MASKS)
            set(hImgCurrentRGB, 'CData', imgRGB);
            
            set(hImgCurrentControlMask, 'CData', maskControl);
            set(hImgCurrentControlCenter, 'XData', blobControl(2), 'YData', blobControl(1));
            
            set(hImgCurrentPlantMask, 'CData', maskPlant);
            set(hImgCurrentPlantCenter, 'XData', blobPlant(2), 'YData', blobPlant(1));
            
            drawnow();
        end
    end
    
    if ((exist('DEBUG_MASKS', 'var') == 1) && DEBUG_MASKS)
        close(hFigCurrentFrame);
    end
end
% --------------------Functions--------------------------------------

function [hueControl, huePlant] = promptForMaskCalibrationValues(imgRGB)
    imgHSV = rgb2hsv(imgRGB);
    
    hFigCalibration = figure();
    axFigCalibrationRGB = subplot(1, 3, 1, 'Parent', hFigCalibration);
    axFigCalibrationControllerMask = subplot(1, 3, 2, 'Parent', hFigCalibration); hold(axFigCalibrationControllerMask, 'on');
    axFigCalibrationPlantMask = subplot(1, 3, 3, 'Parent', hFigCalibration); hold(axFigCalibrationPlantMask, 'on');
    
    imshow(imgRGB, 'Parent', axFigCalibrationRGB);
    
    mask_answer = 'No';
    while(strcmpi(mask_answer, 'No'))
        cla(axFigCalibrationControllerMask);
        cla(axFigCalibrationPlantMask);
        set(axFigCalibrationControllerMask, 'Visible', 'off');
        set(axFigCalibrationPlantMask, 'Visible', 'off');
        
        title(axFigCalibrationRGB, 'Click Control Color');
        [cControl, rControl] = ginput(1);  cControl = round(cControl); rControl = round(rControl);
        rgbControl = double(squeeze(imgRGB(rControl, cControl, :)))' / 255;
        hsvControl = rgb2hsv(rgbControl);
        
        title(axFigCalibrationRGB, 'Click Plant Color');
        [cPlant, rPlant] = ginput(1); cPlant = round(cPlant); rPlant = round(rPlant);
        rgbPlant = double(squeeze(imgRGB(rPlant, cPlant, :)))' / 255;
        hsvPlant = rgb2hsv(rgbPlant);
        
        hueMaskControl = getHueMaskHSV(imgHSV, hsvControl(1));
        hueMaskControl = FillMaskHoles(hueMaskControl);
        hueBlobControlRC = GetBiggestMaskBlob(hueMaskControl);
        
        hueMaskPlant = getHueMaskHSV(imgHSV, hsvPlant(1));
        hueMaskPlant = FillMaskHoles(hueMaskPlant);
        hueBlobPlantRC = GetBiggestMaskBlob(hueMaskPlant);
        
        imshow(hueMaskControl, 'Parent', axFigCalibrationControllerMask);
        plot(axFigCalibrationControllerMask, hueBlobControlRC(2), hueBlobControlRC(1), 'gp');
        
        imshow(hueMaskPlant, 'Parent', axFigCalibrationPlantMask);
        plot(axFigCalibrationPlantMask, hueBlobPlantRC(2), hueBlobPlantRC(1), 'gp');
        
        set(axFigCalibrationControllerMask, 'Visible', 'on');
        set(axFigCalibrationPlantMask, 'Visible', 'on');
        
        mask_answer = questdlg('Are you happy with these color masks?', 'Continue?', 'Yes', 'No', 'Yes');
    end
    
    close(hFigCalibration);
    
    hueControl = hsvControl(1);
    huePlant = hsvPlant(1);
end

function hueMask = getHueMaskHSV(imgHSV, HThresh, HDist)
    if (nargin < 3)
        HDist = 0.02;
    end
    
    SAT_THRESHOLD = 0.5;    % Only pick actual saturated colors
    VAL_THRESHOLD = 0.2;
    
    imgH = imgHSV(:,:,1);
    HLow = HThresh - HDist;
    HHigh = HThresh + HDist;
    
    % Start with any colours that meet our saturation threshold
    satMask = imgHSV(:,:,2) >= SAT_THRESHOLD;
    valMask = imgHSV(:,:,3) >= VAL_THRESHOLD;
    
    hueMask = true(size(imgHSV, 1), size(imgHSV, 2));
    if ((HLow >= 0) && (HHigh <= 1))
        % We can do a simple threshold with a single band
        hueMask = hueMask & ((imgH >= HLow) & (imgH <= HHigh));
    else
        % We have to wrap around to the other side of the Hue channel
        
        if (HLow < 0)
            HLowWrap = mod(HLow, 1);
            hueMaskLow = ((imgH >= HLowWrap) & (imgH <= HThresh));
        else
            hueMaskLow = (imgH >= HLow);
        end
        
        if (HHigh > 1)
            HHighWrap = mod(HHigh, 1);
            hueMaskHigh = ((imgH <= HHighWrap) | (imgH >= HThresh));
        else
            hueMaskHigh = (imgH <= HHigh);
        end
        
        hueMask = hueMaskLow | hueMaskHigh;
    end
    
    hueMask = hueMask & satMask & valMask;
end

function maskOut = FillMaskHoles(maskIn)
    SMALLEST_AREA_THRESHOLD = 100; % Keep areas only if they're bigger than this.
    maskOut = bwareaopen(maskIn, SMALLEST_AREA_THRESHOLD);
    
    % Smooth the border using a morphological closing operation, imclose().
    structuringElement = strel('disk', 4);
    maskOut = imclose(maskOut, structuringElement);
    
    % Fill in any holes in the regions, since they are most likely red also.
    maskOut = imfill(maskOut, 'holes');
end

function blobRC = GetBiggestMaskBlob(maskIn)
    CC = bwconncomp(maskIn);
    
    if (CC.NumObjects == 1)
        idx = 1;
    elseif (CC.NumObjects == 0)
        % There are no blobs in this frame.
        warning('No blob detected');
        blobRC = [];
        return;
    else
        [~, idx] = max(cellfun(@numel, CC.PixelIdxList));
    end
    
    [rows,cols]=ind2sub(CC.ImageSize, CC.PixelIdxList{idx});
    blobRC = [mean([max(rows) min(rows)]), mean([max(cols) min(cols)])];    % Middle of mounding box
end
