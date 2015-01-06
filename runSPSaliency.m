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
I = imread('Tiger.jpg');
%I = img;
%I = cat(3, my_img, my_img, my_img);
%I = imread('boat.jpg'); I = cat(3, I, I, I);
%I = uint8(randn(512,512,3)*255);

%% detect and display superpixels (see spDetect.m)
[E,~,~,segs]=edgesDetect(I,model);
tic, [S,V] = spDetect(I,E,opts); toc
figure; im(I); figure; im(V);

%% compute saliency using superpixels
salMap = SPSaliency(I,S,1);
% a = 5;

%% scale space superpixels
% scale = 1/sqrt(2);
% size = [3 3];
% for i = 1:5
%     h = fspecial('gaussian', size, scale);
%     I = imfilter(I,h);
%     
%     [E,~,~,segs]=edgesDetect(I,model);
%     tic, [S,V] = spDetect(I,E,opts); toc
%     if i == 5
%     figure; im(I); figure; im(V);
%     end
%     scale = scale * sqrt(2);
%     size = size + 2;
%     
%     if i == 3 || i == 5
%         salMap = SPSaliency(I,S, i+1);
%     end
% end


%% compute saliency using superpixels
%salMap = SPSaliency(I,S);
% a = 5;

