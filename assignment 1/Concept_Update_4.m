tic
% Clear workspace

clear all
close all
for p=1:1
    for weighttest=1:25
% Import training data
[TrainingSamplesOfFeatures,TrainingTargetsOfSamples] = importData();

% Numer of neurons in the layers
input=size(TrainingSamplesOfFeatures,2);
outputNeurons=7;
hiddenNeurons=23;
Neurons = [input hiddenNeurons outputNeurons];

% Set learning rate, threshold and maximum number of iterations
alpha=0.8;
beta=0.95;
epochs = 40;

% Trainlength
trainL=round(size(TrainingSamplesOfFeatures,1)*(2/3));

% Creating matrices for efficiency
mse=zeros(trainL,1);
e_k=zeros(1,outputNeurons);

    for i=1:length(Neurons)
        y{i}=zeros(1,Neurons(i));
%     delta{i-1}=zeros(1,Neurons(i));
    end

    for i=2:length(Neurons)
    % Generate random perceptron weights
        w{i-1} = rand(Neurons(i-1),Neurons(i)) -1 * 2;

    % Thresholds
        theta{i-1}= zeros(Neurons(i),1);
    end
% Define sigmoid activation function
sigmoid = @(x)(1/(1+exp(-x)));

% Storage for epoch errors
errors = zeros(epochs, 1);



% Repeat for each epoch
    for epoch = 1:epochs

    
    
    % Alpha gets adjusted with changes in the error
        if epoch>2
        
        % if the error change is small nothing happens if the error change
%         is big alpha increases of the error increases alpha reduces.
        
            if  errors(epoch-2)/errors(epoch-1)<1 || errors(epoch-2)/errors(epoch-1)>1.05
      
%                 alpha=alpha*((errors(epoch-2)/errors(epoch-1)-0.1*alpha));

            end
        
        % if alpha becomes too big reduce it to 2.5
            if alpha>2.5
%             alpha =2.5;
            end
        
%         if alpha<0.05
%             alpha=0.2;
%         end
        
%     else
%         alpha = alphastart;
        end
        
    
    % All combinations of inputs in an epoch
        for inputs = 1:trainL
        
        
        
        % Determine the input
            y{1} = TrainingSamplesOfFeatures(inputs, :);
        
        % Determine the desired output
            y_d = zeros(7,1);
            y_d(TrainingTargetsOfSamples(inputs),1) = 1;
        
        % Determine the output of the hidden layer
        
            for i=1:length(w)
        
                for j = 1:Neurons(i+1)
                
                    y{i+1}(j)=sigmoid(sum(y{i}'.*w{i}(:,j))-theta{i}(j));
      
                end
            end
    
        
        % Determine the error gradient of the output layer

            e_k = y_d - y{length(y)}';
            delta{length(w)} = y{length(y)}.* (1 - y{length(y)}) .* e_k';
            
        % Set new weights to the neurons
            for i=length(w):-1:1
                if i~=length(w)
                    for j = 1:Neurons(i+1)
            
                        delta{i}(j)=y{i+1}(j)*(1-y{i+1}(j)) * sum(delta{i+1}.*w{i+1}(j,:));
                    end
            end
        
                for j = 1:Neurons(i)
                    for k = 1:Neurons(i+1)
                        deltaw(inputs,i,j,k)=alpha * y{i}(j) * delta{i}(k);
                        if inputs>1
                            w{i}(j,k) = w{i}(j,k) + deltaw(inputs,i,j,k)+beta*deltaw(inputs-1,i,j,k);
                        else
                            w{i}(j,k) = w{i}(j,k) + deltaw(inputs,i,j,k);
                        end
                    end
                end
            end
       
        % Calculate the mean squared error
  
            mse(inputs) = sum((1/Neurons(end))*(y_d - y{length(y)}').^2);
  
              
        end
            errors(epoch) = sum(mse);
            disp(['epoch : ' num2str(epoch) ' | error : ' num2str(errors(epoch)) ' | alpha : ' num2str(alpha)])
    % If the total error is zero
        if errors(epoch) == 0
        % Break out the for loop and stop learning
            break
        end
    
        Test=zeros((i-trainL),2);

        clear Test desiredtest actualtest outtest

        if trainL~=length(TrainingSamplesOfFeatures)
    
            desiredTest=TrainingTargetsOfSamples(trainL+1:length(TrainingTargetsOfSamples));    
            actualtest{1}=TrainingSamplesOfFeatures(trainL+1:end,:);
            for i=trainL+1:length(TrainingSamplesOfFeatures)
            
        
                for k=1:length(w)
                    for j=1:Neurons(k+1)
    
                    actualtest{k+1}(i-trainL,j)=sigmoid(dot(actualtest{k}(i-trainL,:),w{k}(:,j))-theta{k}(j,1));
    
                    end
        
                end
        
                [Test(i-trainL,1),Test(i-trainL,2)]=max(actualtest{end}(i-trainL,:));
                
            end
            
        
        end
        results(epoch,weighttest,p)=(sum(round(Test(:,2))==desiredTest)/(length(TrainingTargetsOfSamples)-trainL))*100;
        disp(['Number of correct targets : ' num2str(sum(round(Test(:,2))==desiredTest)) ' which is ' num2str(results(epoch,weighttest,p)) '% of total tested'])


    	

        
    end
    end
    
    X=rand();
    Y=rand();
    Z=rand();
    
%     plot(errors, 'Color',[X,Y,Z],'LineWidth',5)
    hold on
    xlim([1, epoch+1]);
    ylim([0, max(errors)+0.1]);
    figure(2)
    plot(weighttest,max(max(results(:,:,p))), 'd', 'Color',[X,Y,Z],'LineWidth',2)
    hold on


% Set the plot axis limits


Weights{p}=w;
test(p,1)=p;
test(p,2)=errors(end);
    
end
toc