function [motionHist, motionEdges, motionMeanCart] = populateMotion(frame)
%POPULATEMOTION Populates the motion histogram, the edges used, and the
%quantized mean motion vector which represents each bin's motion
%NB: Radius is always first column, meaning the motion histogram goes
%radius down and theta across

flowCart = frame.flowCart;
flowPolar = frame.flowPolar;
binCount = frame.motionBinCount;

% compute histogram and normalize
[motionHist, motionEdges, ~, ~] = histcn(flowPolar, binCount, binCount);
motionHist = motionHist/sum(motionHist(:));

% compute mean motion vector for each bin using cartesian coords
meanU =  histcn(flowPolar, binCount, binCount, 'AccumData', flowCart(:,1), 'Fun', @mean);
meanV =  histcn(flowPolar, binCount, binCount, 'AccumData', flowCart(:,2), 'Fun', @mean);
motionMeanCart = [meanU(:) meanV(:)];

end

