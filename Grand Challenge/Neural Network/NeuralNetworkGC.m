clear all

load('Neurons');

% Define sigmoid activation function
sigmoid = @(x)(1/(1+exp(-x)));

[Samples] = importDataGC();


%% Class for determining the class of a product

for l = 1:2
    
    load(['BestWeights' num2str(l)]);

    actual{1}=Samples(1:length(Samples),:);
        for i=1:length(Samples)
 
            % looping through all the layers and neurons
            for k=1:length(bestSystem)
                for j=1:Neurons(k+1)
                    % calculating the outcomes of the neuron in
                    % specific layer k
                    actual{k+1}(i,j)=sigmoid(dot(actual{k}(i,:),bestSystem{k}(:,j)));

                end

            end     

        end
        
        [f(l,:),g(l,:)]=max(actual{end}');
end