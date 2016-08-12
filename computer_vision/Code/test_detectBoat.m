clear all
close all
%% Read Image
image_o = imread('camera01/MFDC6965.JPG'); % original image
img_g = rgb2gray(image_o);
% figure('Name','Original Image');
% imshow(image_o);
%% kMeans Cluster By Illumination And Colors
image_r = remove(image_o);
cform = makecform('srgb2lab');
lab_img = applycform(image_r, cform);
i_thresVal = 0.05; % threshold to suppress illumination weight
lab_img(:,:,1) = i_thresVal*lab_img(:,:,1);
ab = double(lab_img(:,:,1:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab, nrows*ncols, 3);
nColors = 3; % number of cluster centroids for kMeans
[cluster_idx, cluster_center] = kmeans(ab, nColors, ...
    'distance', 'sqEuclidean', 'Replicates', 6);
pixel_labels = reshape(cluster_idx, nrows, ncols);
% figure('Name', 'Image labeled by cluster index');
% imshow(pixel_labels,[]);
segmented_images = cell(1:nColors);
rgb_label = repmat(pixel_labels, [1 1 nColors]);
for k = 1:nColors
    color = image_r;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end
figure('Name','Clusters by k-Means');
subplot(1,3,1);
imshow(segmented_images{1});
subplot(1,3,2);
imshow(segmented_images{2});
subplot(1,3,3);
imshow(segmented_images{3});

%% Draw Boat Candidate With Boundaries
boat_index = selectBoatImage(segmented_images, nColors);
img_f = rgb2gray(segmented_images{boat_index});
figure('Name','Image Filtered with Boundaries');
imshow(img_f);
hold on
B_f = bwboundaries(img_f); % find boundaries
[rB, cB] = size(B_f);
B_f_size = zeros([rB 1]);
for i = 1:rB
    B_f_size(i) = size(B_f{i}, 1);
end
[sizeVal, B_f_index] = sort(B_f_size, 'descend');
objNum = 9; % number of significant objects
objStartMat = [];
objPosMat = [];
objSizeMat = [];
for k = 1:objNum
plot(B_f{B_f_index(k)}(:,2), B_f{B_f_index(k)}(:,1), 'r-');
hold on
bxL = min(B_f{B_f_index(k)}(:,2)); % left boundary
bxR = max(B_f{B_f_index(k)}(:,2)); % right boundary
byT = min(B_f{B_f_index(k)}(:,1)); % top boundary
byB = max(B_f{B_f_index(k)}(:,1)); % bottom boundary
img_pos = [round((bxL+bxR)/2), round((byT+byB)/2)]; % object position in image coordinates
img_size = [bxR-bxL, byB-byT]; % size measured using boundary rectangle
objStartMat = cat(1, objStartMat, [bxL byT]);
objPosMat = cat(1, objPosMat, img_pos);
objSizeMat = cat(1, objSizeMat, img_size);
rectangle('Position',[bxL byT img_size], 'EdgeColor','y');
end
hold off
%% Prepare Boat Candidates For Recognization
figure('Name','Significant Objects In Original Image');
imshow(image_o);
hold on
for k = 1:objNum
    rectangle('Position',[objStartMat(k,:) objSizeMat(k,:)], 'EdgeColor','y');
    text(objPosMat(k,1), objPosMat(k,2), num2str(k), 'Color','y'); % add labels
    hold on
end
hold off

run(fullfile('build_database_by_ml.m')); 

run(fullfile('new_test.m'));

objScoreMat = [];
objLabelMat = [];
figure('Name', 'Segment Crop');
for k = 1:objNum
    start_r = objStartMat(k, 2);
    start_c = objStartMat(k, 1);
    end_r = start_r + objSizeMat(k, 2);
    end_c = start_c + objSizeMat(k, 1);
    img_obj = image_o(start_r:end_r, start_c:end_c, :);
    [labelIdx, score] = predict(categoryClassifier, img_obj);
    objLabelMat = cat(1, objLabelMat, labelIdx);
    objScoreMat = cat(1, objScoreMat, score);
    subplot(3,3,k);
    imshow(img_obj);
end

[maxSC, maxId] = max(objScoreMat(:,1));
maxId


