% Clear workspace
clear all;

% Set learning rate, threshold and maximum number of iterations
alpha = 0.1;
threshold = 1;
epochs = 50;
hiddenNeurons=12;
outNeurons=7;

% Generate random perceptron weights
weights1=rand(10,12) -1 *2;
weights2=rand(12,7) -1*2;


% Thresholds
thresholds1=rand(12,1) -1 *2;
thresholds2=rand(7,1) -1*2;

% Import training data
[TrainingSamplesOfFeatures,TrainingTargetsOfSamples] = importData();
TrainingSamplesOfFeatures=TrainingSamplesOfFeatures';


% Define sigmoid activation function
sigmoid = @(x)(1/(1+exp(-x)));

% Storage for epoch errors
errors = zeros(epochs, 1);

% Repeat for each epoch
for epoch = 1:epochs
    % All combinations of inputs
    for inputs = 1:size(TrainingSamplesOfFeatures,1)
        
        % Determine input and desired output
        input = TrainingSamplesOfFeatures(:,inputs);
        desired = zeros(1,7);
        desired(1,TrainingTargetsOfSamples(inputs))=1;
        
        for j = 1:hiddenNeurons
            
            % Outputs of hidden neurons           
            actual(j)=sigmoid(dot(input,weights1(:,j))-thresholds1(j,1));
            
        end
            
        for k = 1:outNeurons
            
            % Outputs of output neurons
            out(k)=sigmoid(dot(actual,weights2(:,k))-thresholds2(k,1));
            
        end
           
        % calculating error and delta K
        error = sum(desired-out);
        deltaK = out.*(1-out).*error;
        
        % calculating new weights output layer
        for j=1:hiddenNeurons
            for k=1:outNeurons
                
                weights2(j,k)=weights2(j,k) + alpha*actual(j)*deltaK(k);
                
            end
        end
        
        % Calculating new weights output layer
        deltaJpart=0;
        
        for k=1:outNeurons
            deltaJpart=deltaJpart+deltaK(k).*weights2(:,k);
        end
        
        deltaJ=actual.*(1-actual).*deltaJpart';
        
        for i=1:length(input)
            for j=1:hiddenNeurons
                
                weights1(i,j)=weights1(i,j) + alpha*input(i)*deltaJ(j);
                
            end
        end
                
        
%         % Determine input and desired output
%         input = TrainingSamplesOfFeatures(j,:);
%         desired = TrainingTargetsOfSamples(j);
%         
%         % Calculate the actual perceptron output
%         actual = step(dot(input, weights) - threshold);
%         
%         % Calculate the error
%         error = desired - actual;
%         
%         % Update weights according to the delta rule
%         weights = weights + alpha * input * error;
%         
%         % Add error to total error
%         errors(i) = errors(i) + error * error;
        
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

