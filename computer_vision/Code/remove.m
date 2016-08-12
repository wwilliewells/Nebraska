function img = remove(image)
% remove the unrelated content in the image
[ro, co, zo] = size(image);
img = zeros([ro, co, zo]);
img = im2uint8(img);
btCut = round(ro*0.1);
%r_start = round(ro/4);
r_start = round(1);
r_end = round(ro-btCut);
for i = 1:zo
    img(r_start:r_end,:,i) = image(r_start:r_end,:,i);
end
