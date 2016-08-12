filename = 'image_000';
type = '.jpg';
for i = 1:39
    s = num2str(i);
    fullname = strcat(filename,s,type);
    I = imread(fullname);
    B = imresize(I, [230 370]);
    imwrite(B,fullname);
end
