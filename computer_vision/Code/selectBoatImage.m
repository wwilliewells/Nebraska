function index = selectBoatImage(SegImgMat, SegNum)
% find the index of image containing boat
meanValMat = zeros([SegNum, 1]);
for i = 1:SegNum
    image = SegImgMat{i};
    meanValMat(i) = mean(mean(mean(image(:,:,:))));
end
[minVal, index] = min(meanValMat);