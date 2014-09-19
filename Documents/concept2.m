% Clear workspace
clear all;

% Set learning rate, threshold and maximum number of iterations
alpha = 0.1;
epochs = 50;

% Numer of neurons in the layers
hiddenNeurons = 12;
outputNeurons = 7;

% Generate random perceptron weights
w_ij = rand(10,12) -1 * 2;
w_jk = rand(12,7) - 1 * 2;

% Thresholds
theta_j = ones(12,1) * - 1;
theta_k = ones(7,1) * - 1;

% Import training data
[TrainingSamplesOfFeatures,TrainingTargetsOfSamples] = importData();

% Define sigmoid activation function
sigmoid = @(x)(1/(1+exp(-x)));

% Storage for epoch errors
errors = zeros(epochs, 1);

% Repeat for each epoch
for epoch = 1:epochs
    % All combinations of inputs in an epoch
    for inputs = 1:size(TrainingSamplesOfFeatures,1)
        
        % Determine the input
        x_i = TrainingSamplesOfFeatures(inputs, :);
        
        % Determine the desired output
        y_d_k = zeros(7,1);
        y_d_k(TrainingTargetsOfSamples(inputs),1) = 1;
        
        % Determine the output of the hidden layer
        for j = 1:hiddenNeurons
            
            dotprod = 0;
            for i = 1:size(x_i,2)
                dotprod = dotprod + x_i(i) * w_ij(i,j);
            end
            y_j(j) = sigmoid(dotprod - theta_j(j));
      
        end
        
        % Determine the output of the output layer
        for k = 1:outputNeurons
            
            dotprod = 0;
            for j = 1:size(y_j,2)
                dotprod = dotprod + y_j(j) * w_jk(j,k);
            end
            y_k(k) = sigmoid(dotprod - theta_k(k));
            
        end
        
        % Determine the error gradient of the output layer
        for k = 1:outputNeurons
            
            e_k(k) = y_d_k(k) - y_k(k);
            delta_k(k) = y_k(k) * (1 - y_k(k)) * e_k(k);
            
        end
        
        % Set new weights of the output neurons
        for j = 1:hiddenNeurons
            for k = 1:outputNeurons
                w_jk(j,k) = w_jk(j,k) + alpha * y_j(j) * delta_k(k);
            end
        end
        
        % Determine error gradient of the neurons in the hidden layer
        for j = 1:hiddenNeurons
            
            dotprod = 0;
            for k = 1:size(y_k,2)
                dotprod = dotprod + delta_k(k) * w_jk(j,k);
            end
            delta_j(j) = y_j(j) * (1 - y_j(j)) * dotprod;
            
        end
        
        % Set new weights in the hidden layer
        for i = 1:size(x_i,2)
            for j = 1:hiddenNeurons
                w_ij(i,j) = w_ij(i,j) + alpha * x_i(i) * delta_j(j);
            end
        end
        
        % Calculate the mean squared error
        sum = 0;
        for k = 1:outputNeurons
            sum = sum + (y_d_k(k) - y_k(k))^2;
        end
        errors(epoch) = errors(epoch) + sum / outputNeurons;
              
    end

    % If the total error is zero
    if errors(epoch) == 0
        % Break out the for loop and stop learning
        break
    end
    
end

% Plot the range of errors until epoch i
plot(errors, 'bd-');

% Set the plot axis limits
xlim([1, epoch+1]);
ylim([0, max(errors)+0.1]);

