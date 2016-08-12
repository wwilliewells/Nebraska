%Load Image Sets
imgSets = [ imageSet(fullfile('trainImageSets', 'bird')), ...
            imageSet(fullfile('trainImageSets', 'boat'))];

% display all labels on one line
imgSets.Description
% show the corresponding count of images
imgSets.Count        

%Prepare Training and Validation Image Sets
% determine the smallest amount of images in a category
minSetCount = min([imgSets.Count]);

% Use partition method to trim the set.
imgSets = partition(imgSets, minSetCount, 'randomize');

% Notice that each set now has exactly the same number of images.
imgSets.Count

[trainingSets, validationSets] = partition(imgSets, 0.3, 'randomize');

bird = read(trainingSets(1),1);
boat     = read(trainingSets(2),1);

figure

subplot(1,2,1);
imshow(bird)
subplot(1,2,2);
imshow(boat)


%Create a Visual Vocabulary and Train an Image Category Classifier
bag = bagOfFeatures(trainingSets);

img = read(imgSets(1), 1);
featureVector = encode(bag, img);

% Plot the histogram of visual word occurrences
figure
bar(featureVector)
title('Visual word occurrences')
xlabel('Visual word index')
ylabel('Frequency of occurrence')

categoryClassifier = trainImageCategoryClassifier(trainingSets, bag);


%Evaluate Classifier Performance
confMatrix = evaluate(categoryClassifier, trainingSets);

confMatrix = evaluate(categoryClassifier, validationSets);

% Compute average accuracy
mean(diag(confMatrix));






