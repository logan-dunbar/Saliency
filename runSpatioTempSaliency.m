function temporalSaliencyMap = runSpatioTempSaliency

beginTic = tic;

%% input
if ~exist('input', 'var')
    input.numFrames = 3;
    obj = VideoReader('D:\Logan\Documents\MATLAB\MastersProject\Videos\crane_destroying_building.mp4');
    largeVid = read(obj, [1 input.numFrames]);
    input.vid = uint8(Helper.ResizeVideo(largeVid, 0.5)); % make it slightly more manageable
    input.vidSz = size(input.vid);
end


%% load pre-trained edge detection model and set opts (see edgesDemo.m)
if ~exist('model','var')
    model=load('models/forest/modelBsds');
    model=model.model;
    model.opts.nms=-1;
    model.opts.nThreads=4;
    model.opts.multiscale=0;
    model.opts.sharpen=2;
end

%% set up opts for spDetect (see spDetect.m)
if ~exist('opts','var')
    opts = spDetect;
    opts.nThreads = 2;  % number of computation threads
    opts.k = 2048; %input.vidSz(1) * input.vidSz(2) / 200; % controls scale of superpixels (big k -> big sp) %TODO: tweak k
    opts.alpha = .5;    % relative importance of regularity versus data terms
    opts.beta = .9;     % relative importance of edge versus color terms
    opts.merge = 0;     % set to small value to merge nearby superpixels at end
    %opts.bounds = 0;
end

disp([datestr(datenum(0,0,0,0,0,toc(beginTic)),'HH:MM:SS.FFF') ' <- Input + Model Gen']);

%% compute saliency
framePrev = [];
temporalSaliencyMap = zeros(input.vidSz(1), input.vidSz(2), input.numFrames);
for t = 1:input.numFrames
    tic
    [frameNow, temporalSaliencyMap(:,:,t)] = populateFrame(input.vid(:,:,:,t), framePrev, model, opts, t);
    
    
    
    framePrev = frameNow;
    disp([datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS.FFF') ' <- Frame ' num2str(t)]);
end

%TODO: Super unsafe, but is debug code so that makes it ok.
temporalSaliencyMap = temporalSaliencyMap(:,:,2:end);

disp([datestr(datenum(0,0,0,0,0,toc(beginTic)),'HH:MM:SS.FFF') ' <- Total']);

% %% compute saliency
% im1 = imresize(vid(:,:,:,1),0.5);
% [E1,~,~,~]=edgesDetect(im1,model);
% [S1,~] = spDetect(im1,E1,opts);
% opts.seed = S1; % use superpixels to seed next frame
% 
% for t = 2:nFrames
%     im2 = imresize(vid(:,:,:,t),0.5);
%     [E2,~,~,~]=edgesDetect(im2,model);
%     [S2,~] = spDetect(im2,E2,opts);
%     opts.seed = S2;
%     
%     % color histograms (CH{1} is global)
%     im2Lab = RGB2Lab(im2);
%     im2Lab = reshape(im2Lab, [], 3);
%     
%     sps.inds = unique(S2);
%     sps.inds = sps.inds(sps.inds > 0);
%     sps.n = length(sps.inds);
%     
%     CH = cell(sps.n + 1, 1);
%     [CH{1}.hist, ~, CH{1}.mid, CH{1}.loc] = histcn(im2Lab, l_bin_edges, ab_bin_edges, ab_bin_edges);
%     
%     for i = 2:sps.n+1
%         sps.loc = S2 == sps.inds(i-1);
%         sps.curInds = find(sps.loc);
%         sps.vals = im2Lab(sps.curInds,:);
%         [y,x] = ind2sub(imSz, inds);
%         CH{i}.center = 
%         
%         [CH{i}.hist, ~, CH{i}.mid, CH{i}.loc] = histcn(sps.vals, l_bin_edges, ab_bin_edges, ab_bin_edges);
%         
%         
%         
%     end
%     
%     
%     
%     
%     im1 = im2;
%     S1 = S2;
% end
% 
% end
% 
% function [S, colorHist, colorMid, colorLoc] = computeColorHistogram(imgIn, q, model, opts)
% 
% 
% 
% % convert to L*a*b*
% labImg = RGB2Lab(imgIn);
% labImg = reshape(labImg, [], 3);
% 
% % quantisation edges
% l_bin_edges = 0:100/q:100;
% ab_bin_edges = -128:255/q:127;
% 
% spsInds = unique(spsMap); % unique and sort (built in)
% spsInds = spsInds(spsInds > 0);
% spsTotCount = length(spsInds);
% % for i
% 
% 
% % compute global color histogram
% [colorHist, ~, colorMid, colorLoc] = histcn(labImg, l_bin_edges, ab_bin_edges, ab_bin_edges);
% colorHist = colorHist/sum(colorHist(:));
% 
% end