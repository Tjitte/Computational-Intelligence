
% Clear workspace
clear all
close all

%% initializing
% Set learning rate,momentum rate, number of epochs per test
alpha=1.5;
beta=0.95;
epochs = 20;
errors = zeros(epochs, 1);          % Matrix for total epoch errors
weighttests=100;

% Import training data
[TrainingSamplesOfFeatures,TrainingTargetsOfSamples,SamplesOfFeatures] = importData();
answers2=zeros((1/6)*length(TrainingSamplesOfFeatures),7,epochs,weighttests,1);

for p=1:1
    
    for weighttest=1:weighttests

        % Number of neurons in the layers
        input=size(TrainingSamplesOfFeatures,2);
        outputNeurons=7;
        hiddenNeurons=15;
        Neurons = [input hiddenNeurons outputNeurons];

 

        % Trainlength
        trainL=round(size(TrainingSamplesOfFeatures,1)*(2/3));

        % Creating matrices for efficiency
        mse=zeros(trainL,1);                % Matrix for mean squared error
        e_k=zeros(1,outputNeurons);         % Matrix for current error
        
        desiredTest2= zeros(7,1);
        
        for i=2:length(Neurons)
            % Generate random perceptron weights
                w{i-1} = rand(Neurons(i-1),Neurons(i)) -1 * 2;

            % Thresholds
                theta{i-1}= zeros(Neurons(i),1);
        end
        
        % Define sigmoid activation function
        sigmoid = @(x)(1/(1+exp(-x)));

%% Training

        % Repeat for each epoch
        for epoch = 1:epochs



            % Alpha gets adjusted with changes in the error
            if epoch>2

                % if the error change is small alpha gets a little smaller if the error change
                % is big alpha increases if the error change is less than 0
                % alpha gets a lot smaller

                if  errors(epoch-2)/errors(epoch-1)<1 || errors(epoch-2)/errors(epoch-1)>1.05

                    alpha=alpha*((errors(epoch-2)/errors(epoch-1)-0.2*alpha));
                    
                else
                    alpha=alpha*0.98;
                end

                % if alpha becomes too big reduce it to 1
                if alpha>1.5
                    alpha =1.5;
                end

            end


            % for loop loops through all the data in the training set
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

                % looping through all te weights from each layer
                for i=length(w):-1:1
                    % when in the hidden layer->
                    if i~=length(w)
                        % looping through all the neurons in the hidden layer
                        for j = 1:Neurons(i+1)
                            % calculating error gradient with following formula
                            delta{i}(j)=y{i+1}(j)*(1-y{i+1}(j)) * sum(delta{i+1}.*w{i+1}(j,:));

                        end
                    end

                    % looping throught each weight
                    for j = 1:Neurons(i)

                        for k = 1:Neurons(i+1)
                            % creating matrix with delta weights
                            deltaw(inputs,i,j,k)=alpha * y{i}(j) * delta{i}(k);
                            if inputs>1
                                % updating all the weights
                                w{i}(j,k) = w{i}(j,k) + deltaw(inputs,i,j,k)+beta*deltaw(inputs-1,i,j,k);
                            end
                        end
                    end
                end

                % Calculate the mean squared error

                mse(inputs) = sum((e_k).^2);


            end
            
            % sum of all the errors in 1 epoch and weigthtest
            errors(epoch,weighttest,p) = sum(mse)/(trainL*outputNeurons);
            % display information about this epoch
            disp(['epoch : ' num2str(epoch) ' | error : ' num2str(errors(epoch,weighttest,p)) ' | alpha : ' num2str(alpha)])

            % if the training length is shorter than the total number of
            % elements the following happens
            
            %% Validation
            
            if trainL~=length(TrainingSamplesOfFeatures)
                
                % determining the desired outcomes for the validation set
                desiredTest=TrainingTargetsOfSamples(trainL+1:trainL+1+length(TrainingSamplesOfFeatures)/6);
                % inputvalues for system to be tested
                actualtest{1}=TrainingSamplesOfFeatures(trainL+1:end,:);
                
                % looping through all the elements in the validation set
                for i=trainL+1:trainL+1+length(TrainingSamplesOfFeatures)/6
                    % determining the desired outcome
                    desiredTest2(TrainingTargetsOfSamples(i),i-trainL) = 1;
                    % looping through all the layers and neurons
                    for k=1:length(w)
                        for j=1:Neurons(k+1)
                            % calculating the outcomes of the neuron in
                            % specific layer k
                            actualtest{k+1}(i-trainL,j)=sigmoid(dot(actualtest{k}(i-trainL,:),w{k}(:,j))-theta{k}(j,1));
    
                        end
        
                    end
                    
                    % determing of the output layer which element is
                    % biggest and at which location it was
                    [Test(i-trainL,1),Test(i-trainL,2)]=max(actualtest{end}(i-trainL,:));
                    % validating error between desired output and actual
                    % output
                    val(i-trainL,:)=(actualtest{end}(i-trainL,:)-desiredTest2(:,i-trainL)').^2;
                    
                end
                
            end
        
        % results for how much correct guesses
        results(epoch,weighttest,p)=(sum(round(Test(:,2))==desiredTest)/(length(TrainingTargetsOfSamples)/6))*100;
        % results for how big the total summed error was between desired
        % outputs and actual outputs
        results2(epoch,weighttest,p)=sum(sum(val))/(length(TrainingSamplesOfFeatures)/6);
        
        % display information about how good the validation was
        disp(['test: ' num2str(weighttest) ' number of hidden neurons: ' num2str(hiddenNeurons)])
        disp(['Number of correct targets : ' num2str(sum(round(Test(:,2))==desiredTest)) ' which is ' num2str(results(epoch,weighttest,p)) '% of total tested'])
        
        % saving all the weights for each epoch
        Weights{p,weighttest,epoch}=w;
            
        end
        
        % for the neuron set picking the best system accoring to the validation
        [a,b]=min(results2);
        [c,d]=min(a);
        bestSystem=Weights{1,d,b};
        
        save(['BestWeights' num2str(weighttest)] ,'bestSystem');
        
    end

    save('Neurons','Neurons');


    %% Testing
    

    plot(errors(:,d,1),'b')
    hold on
    plot(results2(:,d,1),'r')
    
    % creating a matrix for the output
    d_out=zeros(7,length(TrainingSamplesOfFeatures)/6);
    
    % determining the desired output
    desiredout=TrainingTargetsOfSamples((trainL+length(TrainingSamplesOfFeatures)/6+1):length(TrainingTargetsOfSamples));
    
    % setting actual{1} as the input data
    actual{1}=TrainingSamplesOfFeatures(trainL+length(TrainingSamplesOfFeatures)/6+1:end,:);
        for i=trainL+length(TrainingSamplesOfFeatures)/6+1:length(TrainingSamplesOfFeatures)
            % determining the desired outcome
            d_out(desiredout(i-(trainL+length(TrainingSamplesOfFeatures)/6)),i)=1;
            % looping through all the layers and neurons
            for k=1:length(bestSystem)
                for j=1:Neurons(k+1)
                    % calculating the outcomes of the neuron in
                    % specific layer k
                    actual{k+1}(i-(trainL+length(TrainingSamplesOfFeatures)/6),j)=sigmoid(dot(actual{k}(i-(trainL+length(TrainingSamplesOfFeatures)/6),:),bestSystem{k}(:,j))-theta{k}(j,1));

                end

            end     
            
        end
    


    % Calculating the error on the test set
    for i=1:length(TrainingSamplesOfFeatures)/6
        err_out(i,:)=(actual{end}(i,:)-d_out(i)).^2;
    end
    
    total_test_error=sum(sum(err_out))/length(TrainingSamplesOfFeatures)/6;
    
    % determining which the output for the test set and comparing to the
    % real output data
    [e,f]=max(actual{end}');
    
    answers=zeros(length(TrainingSamplesOfFeatures)/6,7);
    
    for i=1:length(TrainingSamplesOfFeatures)/6
        answers(i,f(i))=1;
    end
    f=f';
    vergelijking=sum(f==desiredout)/(length(TrainingSamplesOfFeatures)/6)*100;
    
    disp(['The error on the test set is: ' num2str(total_test_error) ' | percentage right guessed: ' num2str(vergelijking) ' | best valdidation set: ' num2str(results(b(d),d,1))])
    

end

%% Using the system weights for determining the unknown targets

    actual{1}=SamplesOfFeatures(1:length(SamplesOfFeatures),:);
        for i=1:length(SamplesOfFeatures)
 
            % looping through all the layers and neurons
            for k=1:length(bestSystem)
                for j=1:Neurons(k+1)
                    % calculating the outcomes of the neuron in
                    % specific layer k
                    actual{k+1}(i,j)=sigmoid(dot(actual{k}(i,:),bestSystem{k}(:,j))-theta{k}(j,1));

                end

            end     

        end

% Determining output for unknown set
[f,g]=max(actual{end}');


% Writing to txt file
fid=fopen('36_classes.txt','wt');
for i=1:length(SamplesOfFeatures)
    if i < length(SamplesOfFeatures)
        fprintf(fid,[num2str(g(i)) ', ']);
    else
        fprintf(fid,[num2str(g(i))]);   
    end
end

% creating vectors for confusion matrix
targets=zeros(length(TrainingSamplesOfFeatures)/6,7);
targets_test=TrainingTargetsOfSamples(trainL+length(TrainingSamplesOfFeatures)/6+1:length(TrainingSamplesOfFeatures));
for i=1:length(TrainingSamplesOfFeatures)/6
    targets(i,targets_test(i))=1;
end

% making plots look nice
figure(2)
plotconfusion(targets',answers');
figure(1)
xlabel('epoch -->')
ylabel('mean squared error -->')
title('training error and validation error per epoch  (red validation error - blue training error)')
