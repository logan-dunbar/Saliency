function [spArray, motionSignificanceHistory, saliencySpatial] = populateSpArray(frameNow, framePrev)
%POPULATESPARRAY Populates the superpixel array by looping through each superpixel and computing its
%histograms and properties

% motion significance window length
motionSignigicanceWinLength = 5;
motionSignificanceHistory = [];
saliencySpatial = zeros(frameNow.imgSize);

% preallocate spArray
spArray = repmat(struct( ...
    'pixelCount',[], ...
    'yInds', [], ...
    'xInds', [], ...
    'center',[], ...
    'pixelVals',[], ...
    'flowVals',[], ...
    'motionHist',[], ...
    'colorHist',[], ...
    'meanMotion',[], ...
    'interFrameSimilarity', [], ...
    'saliencyMotionDistinctiveness', [], ...
    'saliencyTemporalPrediction', [], ...
    'saliencyTemporal', [], ...
    'intraFrameSimilarity', zeros(frameNow.spNum, 1), ...
    'saliencyGlobalContrast', [], ...
    'saliencySpatialSparcity', [], ...
    'saliencySpatial', []), frameNow.spNum, 1);

if frameNow.frameNum ~= 1
    prevCenters = reshape([framePrev.spArray.center], 2, [])';
end

% loop through each superpixel
for i = 1:frameNow.spNum
    % get indices and count
    spArray(i).pixelInds = find(frameNow.spMap == frameNow.spInds(i));
    spArray(i).pixelCount = length(spArray(i).pixelInds);
    
    % compute center
    [spArray(i).yInds, spArray(i).xInds] = ind2sub(frameNow.imgSize, spArray(i).pixelInds);
    spArray(i).center = [sum(spArray(i).yInds)/length(spArray(i).pixelInds) sum(spArray(i).xInds)/length(spArray(i).pixelInds)];
    spArray(i).distanceToFrameCenter = norm(spArray(i).center - frameNow.center);
    
    % compute superpixel level color histogram and normalize
    [spArray(i).colorHist, ~, ~, spArray(i).colorLoc] = histcn(frameNow.labImg(spArray(i).pixelInds,:), frameNow.colorEdges{1}, frameNow.colorEdges{2}, frameNow.colorEdges{3});
    spArray(i).colorHist = spArray(i).colorHist/sum(spArray(i).colorHist(:));
    
    % compute superpixel level motion histogram and normalize
    if frameNow.frameNum ~= 1
        spArray(i).motionHist = histcn(frameNow.flowPolar(spArray(i).pixelInds,:), frameNow.motionEdgesPolar{1}, frameNow.motionEdgesPolar{2});
        spArray(i).motionHist = spArray(i).motionHist/sum(spArray(i).motionHist(:));
        
        % calculate mean flow of the superpixel
        spMeanFlow = mean(frameNow.flowCart(spArray(i).pixelInds,:), 1);
        % project the center backwards
        projectedCenter = repmat(spArray(i).center - spMeanFlow, framePrev.spNum, 1);
        % compute distances between projected center and all previous superpixel centers
        norms = sqrt(sum((projectedCenter - prevCenters).^2, 2));
        % find the closest match
        [~, minSpIndex] = min(norms);
        % find neighbouring indices of the most likely superpixel from the previous frame
        spArray(i).prevNeighbouringSpInds = findNeighbouringSuperpixelIndices(minSpIndex, framePrev);
        % calculate interframe similarity between current superpixel and each neighbour
        spArray(i).interFrameSimilarity = zeros(length(spArray(i).prevNeighbouringSpInds), 1);
        currSp = spArray(i);
        for j = 1:length(spArray(i).prevNeighbouringSpInds)
            prevSp = framePrev.spArray(spArray(i).prevNeighbouringSpInds(j));
            spArray(i).interFrameSimilarity(j) = sum(sqrt(currSp.colorHist(:) .* prevSp.colorHist(:))) * ...
                (1 - norm(currSp.center - spMeanFlow - prevSp.center)/framePrev.diagonal);
        end
    end
end

distancesToFrameCenter = [spArray.distanceToFrameCenter]';

for i = 1:frameNow.spNum
    %     spArray(i).saliencyGlobalContrast = 0;
    %     for j = 1:frameNow.contrastBinCount^3
    %         innerSum = 0;
    %         for k = 1:frameNow.contrastBinCount^3
    %             innerSum = innerSum + norm(frameNow.colorMean(j,:) - frameNow.colorMean(k,:)) * frameNow.colorHist(k);
    %         end
    %         spArray(i).saliencyGlobalContrast = spArray(i).saliencyGlobalContrast + (spArray(i).colorHist(j) * innerSum);
    %     end
    
    % calculate global contrast
    % vectorized form of the global contrast saliency calculation. It flattens out the double loop
    % into large q^3 x q^3 matrices, computes the inner sum for all combinations of j and k,
    % multiplies each value by the correspoing global contrast bin count at k, sums that value, and
    % then does the outer j sum and multiply with the superpixel color bin count. Higher memory
    % requirement but much, much faster. ~175x faster!
    color_l1 = repmat(frameNow.colorMean(:,1), 1, frameNow.contrastBinCount^3);
    color_l2 = repmat(frameNow.colorMean(:,1)', frameNow.contrastBinCount^3, 1);
    color_a1 = repmat(frameNow.colorMean(:,2), 1, frameNow.contrastBinCount^3);
    color_a2 = repmat(frameNow.colorMean(:,2)', frameNow.contrastBinCount^3, 1);
    color_b1 = repmat(frameNow.colorMean(:,3), 1, frameNow.contrastBinCount^3);
    color_b2 = repmat(frameNow.colorMean(:,3)', frameNow.contrastBinCount^3, 1);
    
    color_l_diff = color_l1 - color_l2;
    color_a_diff = color_a1 - color_a2;
    color_b_diff = color_b1 - color_b2;
    
    color_norm = sqrt(color_l_diff.^2 + color_a_diff.^2 + color_b_diff.^2);
    globalColorHist = repmat(frameNow.colorHist(:), 1, frameNow.contrastBinCount^3);
    innerSum = sum(color_norm.*globalColorHist, 1);
    spArray(i).saliencyGlobalContrast = sum(spArray(i).colorHist(:)'.*innerSum);
    
    % TODO: see if I can't vectorize this better, just worried about memory constraints:
    % 700x700x4096 * 2 for the first term...
    spArray(i).intraFrameSimilarity = zeros(frameNow.spNum, 1);
    for j = 1:frameNow.spNum
        spArray(i).intraFrameSimilarity(j) = sum(sqrt(spArray(i).colorHist(:) .* spArray(j).colorHist(:))) * ...
            (1 - norm(spArray(i).center - spArray(j).center)/frameNow.diagonal);
    end
    
    spArray(i).colorSpatialDistribution = sum(spArray(i).intraFrameSimilarity.*distancesToFrameCenter)/sum(spArray(i).intraFrameSimilarity);
end

% calculate saliencySpatialSparcity and saliencySpatial
maxColorSpatialDistribution = max([spArray.colorSpatialDistribution]);
minColorSpatialDistribution = min([spArray.colorSpatialDistribution]);
for i = 1:frameNow.spNum
    spArray(i).saliencySpatialSparcity = (maxColorSpatialDistribution - spArray(i).colorSpatialDistribution)/(maxColorSpatialDistribution - minColorSpatialDistribution);
    spArray(i).saliencySpatial = spArray(i).saliencySpatialSparcity * spArray(i).saliencyGlobalContrast;
end

for i = 1:frameNow.spNum
    % find neighbour indices for each superpixel in the current frame
    spArray(i).nowNeighbouringSpInds = findNeighbouringSuperpixelIndices(i, frameNow, spArray);
    % get the bin indices of each pixel in the superpixel
    colorBinInds = sub2ind(repmat(frameNow.contrastBinCount, 1, 3), spArray(i).colorLoc(:,1), spArray(i).colorLoc(:,2), spArray(i).colorLoc(:,3));
    % compute the pixel-level spatial saliency score for each pixel in the superpixel
    numeratorSum = zeros(spArray(i).pixelCount, 1);
    denomenatorSum = zeros(spArray(i).pixelCount, 1);
    for j = 1:length(spArray(i).nowNeighbouringSpInds)
        numeratorSum = numeratorSum + (spArray(j).saliencySpatial * spArray(j).colorHist(colorBinInds));
        denomenatorSum = denomenatorSum + spArray(j).colorHist(colorBinInds);
    end
    saliencySpatial(spArray(i).pixelInds) = numeratorSum./denomenatorSum;
end


if frameNow.frameNum ~= 1
    % calculate motion distinctiveness
    % vectorized form of the motion distinctiveness calculation. It flattens out the double loop
    % into large q^2 x q^2 matrices, computes the inner sum for all combinations of j and k,
    % multiplies each value by the correspoing global motion bin count at k, sums that value, and
    % then does the outer j sum and multiply with the superpixel motion bin count. Higher memory
    % requirement but much, much faster.
    for i = 1:frameNow.spNum
        flow_u1 = repmat(frameNow.motionMeanCart(:,1), 1, frameNow.motionBinCount^2);
        flow_u2 = repmat(frameNow.motionMeanCart(:,1)', frameNow.motionBinCount^2, 1);
        flow_v1 = repmat(frameNow.motionMeanCart(:,2), 1, frameNow.motionBinCount^2);
        flow_v2 = repmat(frameNow.motionMeanCart(:,2)', frameNow.motionBinCount^2, 1);
        
        flow_u_diff = flow_u1 - flow_u2;
        flow_v_diff = flow_v1 - flow_v2;
        
        flow_norm = sqrt(flow_u_diff.^2 + flow_v_diff.^2);
        
        globalMotionHist = repmat(frameNow.motionHist(:), 1, frameNow.motionBinCount^2);
        innerSum = sum(flow_norm.*globalMotionHist, 1);
        outerSum = sum(spArray(i).motionHist(:)'.*innerSum);
        spArray(i).saliencyMotionDistinctiveness = outerSum;
    end
    
    % calculate motion significance
    motionSignificance = sum([spArray.saliencyMotionDistinctiveness].*[spArray.pixelCount]);
    if frameNow.frameNum == 2
        motionSignificanceHistory = motionSignificance;
    else
        % calculate the adjustment weight based on motion significance for each frame
        adjustmentWeight = max(median(framePrev.motionSignificanceHistory)/(median(framePrev.motionSignificanceHistory) + motionSignificance), 0.5);
        motionSignificanceHistory = [framePrev.motionSignificanceHistory motionSignificance];
        motionSignificanceHistory = motionSignificanceHistory(max(length(motionSignificanceHistory)- motionSignigicanceWinLength, 1):end);
    end
    
    % TODO: double loop, see if I can't collapse it/vectorize it
    for i = 1:frameNow.spNum
        if frameNow.frameNum == 2
            % no previous saliency value on first motion frame, just use motion distinctiveness
            spArray(i).saliencyTemporal = spArray(i).saliencyMotionDistinctiveness;
        else
            if sum(spArray(i).interFrameSimilarity) > 0
                % calculate temporal saliency prediction
                prevTemporalSaliency = [framePrev.spArray(spArray(i).prevNeighbouringSpInds).saliencyTemporal]';
                spArray(i).saliencyTemporalPrediction = sum(spArray(i).interFrameSimilarity.*prevTemporalSaliency)/sum(spArray(i).interFrameSimilarity);
                
                % calculate temporal saliency by adjustment
                spArray(i).saliencyTemporal = adjustmentWeight*spArray(i).saliencyTemporalPrediction + ...
                    (1 - adjustmentWeight)*spArray(i).saliencyMotionDistinctiveness;
            else
                spArray(i).saliencyTemporalPrediction = 0;
                spArray(i).saliencyTemporal = spArray(i).saliencyMotionDistinctiveness;
            end
        end
    end
end

end

