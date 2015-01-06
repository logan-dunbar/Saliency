%% load pre-trained edge detection model and set opts (see edgesDemo.m)
if ~exist('model','var')
    model=load('models/forest/modelBsds'); model=model.model;
    model.opts.nms=-1; model.opts.nThreads=4;
    model.opts.multiscale=0; model.opts.sharpen=2;
end

%% set up opts for spDetect (see spDetect.m)
if ~exist('opts','var')
    opts = spDetect;
    opts.nThreads = 2;  % number of computation threads
    opts.k = 512;       % controls scale of superpixels (big k -> big sp)
    opts.alpha = .5;    % relative importance of regularity versus data terms
    opts.beta = .9;     % relative importance of edge versus color terms
    opts.merge = 0;     % set to small value to merge nearby superpixels at end
    opts.bounds = 0;
end

%% input
obj = VideoReader('D:\Logan\Documents\MATLAB\MastersProject\Videos\crane_destroying_building.mp4');
vid = read(obj, [1 2]);
im1 = imresize(vid(:,:,:,1),0.5);
im2 = imresize(vid(:,:,:,2),0.5);

%% detect and display superpixels (see spDetect.m)
[E,~,~,segs]=edgesDetect(im1,model);
tic, [S,V] = spDetect(im1,E,opts); toc
figure; im(im1); figure; im(V);

[E1,~,~,segs1]=edgesDetect(im2,model);
tic, [S1,V1] = spDetect(im2,E1,opts); toc
figure; im(im2); figure; im(V1);

%% compute saliency using superpixels
if ~exist('flow','var'); disp('motion calc');tic;flow = mex_LDOF(double(im1), double(im2));toc; end
sal = MotionSaliency(flow,S);
