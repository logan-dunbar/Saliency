function [ frameNow, temporalSalMap ] = populateFrame(img, framePrev, model, opts, frameNum)
%POPULATEFRAME Populates the frame level variables for calculating saliency

frameNow.frameNum = frameNum;

% number of bins per color channel
frameNow.contrastBinCount = 16;
% number of bins for magnitude and angle of flow field
frameNow.motionBinCount = 16;

frameNow.img = img;
[frameNow.imgSize(1), frameNow.imgSize(2), ~] = size(img);
frameNow.diagonal = sqrt(sum(frameNow.imgSize.^2));
frameNow.center = frameNow.imgSize./2;

% convert to L*a*b* and reshape
frameNow.labImg = reshape(RGB2Lab(double(img)), [], 3);

% compute superpixels
[edges, ~, ~, ~]=edgesDetect(img, model);
[frameNow.spMap, ~] = spDetect(img, edges, opts);
frameNow.spInds = unique(frameNow.spMap);
frameNow.spInds = frameNow.spInds(frameNow.spInds > 0);
frameNow.spNum = length(frameNow.spInds);

% global contrast histogram
[frameNow.colorHist, frameNow.colorEdges, frameNow.colorMean] = populateContrast(frameNow);

% compute optical flow (For Polar coords, chose to put radius in first column)
if frameNow.frameNum ~= 1
    frameNow.flowCart = reshape(mex_LDOF(double(framePrev.img), double(frameNow.img)), [], 2);
    [frameNow.flowPolar(:,2), frameNow.flowPolar(:,1)] = cart2pol(frameNow.flowCart(:,1), frameNow.flowCart(:,2));
    [frameNow.motionHist, frameNow.motionEdgesPolar, frameNow.motionMeanCart] = populateMotion(frameNow);
end

% popluate the array of superpixel objects with their color and motion histograms etc
[frameNow.spArray, frameNow.motionSignificanceHistory, frameNow.saliencySpatial] = populateSpArray(frameNow, framePrev);

% TODO: debugging code to view the saliency map, remove when not needed
temporalSalMap = zeros(frameNow.imgSize);
if frameNow.frameNum ~= 1
    for i=1:frameNow.spNum
        spsLoc = frameNow.spMap == frameNow.spInds(i);
        temporalSalMap(spsLoc) = frameNow.spArray(i).saliencyTemporal;
    end
end

end

