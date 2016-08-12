
%Try the Newly Trained Classifier on Test Images
%Test Boat
img_boat = imread('boat_test_0001.jpg');
[labelIdx_boat] = predict(categoryClassifier, img_boat);

% Display the string label
categoryClassifier.Labels(labelIdx_boat)
position =  [37 37];
value = categoryClassifier.Labels(labelIdx_boat);
RGB_boat = insertText(img_boat,position,value,'AnchorPoint','LeftBottom');
figure
imshow(RGB_boat),title('Result of boat test');


%Test Bird
img_bird = imread(fullfile('bird_test_0001.jpg'));
[labelIdx_bird] = predict(categoryClassifier, img_bird);

%Display the string label
categoryClassifier.Labels(labelIdx_bird)
position =  [37 37];
value = categoryClassifier.Labels(labelIdx_bird);
RGB_bird = insertText(img_bird,position,value,'AnchorPoint','LeftBottom');
figure
imshow(RGB_bird),title('Result of bird test');
