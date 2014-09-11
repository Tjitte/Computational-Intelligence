% Clear workspace
clear all;

% Set learning rate, threshold and maximum number of iterations
alpha = 0.1;
threshold = 1;
epochs = 100;

% Generate random perceptron weights
weights = rand(1,10) * 2 - 1;

% Import training data
[TrainingSamplesOfFeatures,TrainingTargetsOfSamples] = importData();

% Define step activation function (as an anonymous function, neat trick!)
step = @(x)(x >= 0);

% Storage for epoch errors
errors = zeros(epochs, 1);

% Repeat for each epoch
for i = 1:epochs
    % All combinations of inputs
    for j = 1:size(TrainingSamplesOfFeatures,1)
        
        % Determine input and desired output
        input = TrainingSamplesOfFeatures(j,:);
        desired = TrainingTargetsOfSamples(j);
        
        % Calculate the actual perceptron output
        actual = step(dot(input, weights) - threshold);
        
        % Calculate the error
        error = desired - actual;
        
        % Update weights according to the delta rule
        weights = weights + alpha * input * error;
        
        % Add error to total error
        errors(i) = errors(i) + error * error;
        
    end
    
    % If the total error is zero
    if errors(i) == 0
        % Break out the for loop and stop learning
        break
    end
    
end

% Plot the range of errors until epoch i
plot(errors, 'bd-');

% Set the plot axis limits
xlim([1, i+1]);
ylim([0, max(errors)+0.1]);

