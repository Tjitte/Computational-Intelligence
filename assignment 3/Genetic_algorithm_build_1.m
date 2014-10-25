% clear all;
% 
% %list of variables
% iterations = 1;
% 
% %% reading product locations
% fd = fopen('mazes/tsp products.txt');
% itemNumber = textscan(fd,'%f','delimiter',';:,');
% fclose(fd);
% 
% for i=3:3:length(itemNumber{1,1})-1
%     locations(i/3,1)=(itemNumber{1,1}(i));
% end
% 
% for i=4:3:length(itemNumber{1,1})
%     locations((i-1)/3,2)=((itemNumber{1,1}(i)));
% end
% 
% %% Calculating pathlengths between two points in a maze
% 
% % initializing
% nodes = length(locations(:,1));
% connection = (nodes*(nodes-1))/2;
% connections = zeros(1,nodes);
% f=@shortestPath;
% counter = 0;
% 
% % Calculating paths and placing in connections
% for x=1:connection
%     tempnodes = nodes;
%     while(tempnodes > x)
%         connections(x,tempnodes) = f(locations(x,:),locations(tempnodes,:),1,6);
%         connections(tempnodes,x) = connections(x,tempnodes);
% 
%         tempnodes = tempnodes-1;
%         counter = counter + 1;
%         disp(counter)
%     end
% end

%% Creating random chromosomes
amountchromo = 100;
chromo = zeros(amountchromo,nodes,iterations);
random = 1;

%
for i=1:length(chromo(:,1,1))
    for j=1:length(chromo(1,:,1))
        random = round(rand()*18+0.5);
        while sum(random == chromo(i,:,1)) ~= 0
            random = round(rand()*18+0.5);
        end
        chromo(i,j,1) = random;
    end
end

fitness = zeros(amountchromo,iterations);
fitnessRatio = zeros(amountchromo,iterations);


for iter=1:iterations
    %% fitness ratio
    distribution = zeros(amountchromo,iter);

    for i=1:length(chromo(:,1,iter))
        for j=1:length(chromo(1,:,iter))-1
            point(1) = chromo(i,j,iter);
            point(2) = chromo(i,j+1,iter);
            fitness(i,iter) = fitness(i,iter) + connections(point(2),point(1));
        end
    end
    pathlength=min(fitness(:,iter));
    fitness(:,iter) = 1./fitness(:,iter)*100000;
    fitnessRatio(:,iter) = fitness(:,iter)./sum(fitness(:,iter));

   
    % choosing the new population
    distribution(1,iter)=0;
    distribution(2:amountchromo+1,iter) = cumsum(fitnessRatio(:,iter));
    
    amountselections = 2;
    
    for i=1:amountselections:amountchromo
       %% selecting
       
        for k=1:amountchromo/amountselections
            select =zeros(1,amountselections);
            while length(unique(select)) ~= amountselections
                for p=1:amountselections
                    random = rand;
                    for j=1:amountchromo
                        if distribution(j,iter) < random && distribution(j+1,iter) > random
                            select(p) = j;
                            break
                        end
                    end
                end
            end


        end
    
    %% Crossover

        pcrossover = 0.7;
        random2 = rand;

        if random2 < pcrossover
            random = round(rand()*18+0.5);

            tempvect=zeros(1,nodes);

            tempvect(1:random)=chromo(select(1),1:random,iter);
            
            newchromoloop=random;
            for loopchromo=1:nodes
                
                if sum(chromo(select(2),loopchromo,iter) ==  tempvect) == 0;
                    newchromoloop=newchromoloop+1;
                    tempvect(newchromoloop)= chromo(select(2),loopchromo,iter);
                    
                end
            end

        end
        chromo(i,:,iter+1)=tempvect;
        
        if random2 < pcrossover

            tempvect=zeros(1,nodes);

            tempvect(1:random)=chromo(select(2),1:random,iter);
            
            newchromoloop=random;
            for loopchromo=1:nodes
                
                if sum(chromo(select(1),loopchromo,iter) ==  tempvect) == 0;
                    newchromoloop=newchromoloop+1;
                    tempvect(newchromoloop)= chromo(select(1),loopchromo,iter);
                    
                end
            end
            
            chromo(i+1,:,iter+1)=tempvect;

        end
        
        
        if random2 >= pcrossover
            chromo(i,:,iter+1)=chromo(select(1),:,iter);
            chromo(i+1,:,iter+1)=chromo(select(2),:,iter);
            
        end
    end

    
  
    
    %% mutation
    pmutation = 0.005;
    for i=1:amountchromo
        for j=1:nodes
            random = rand;
            if rand < pmutation
                change = chromo(i,j,iter+1);
                random2=round(rand()*18+0.5);
                chromo(i,j,iter+1)=chromo(i,random2,iter+1);
                chromo(i,random2,iter+1)=change;

            end
        end
    end
    
    disp(['pathlength: ' num2str(pathlength)]);
end

