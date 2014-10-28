tic
clear all;
clc

fid = fopen('Pathwayf.mat');
if fid~=-1
    load('Pathwayf.mat','Pathwayf');
end

%% list of variables
iterations = 100000;
runs = 30;
iters = 30;
Reiterate = 20;

% start and end coordinates
fd = fopen('mazes/hard coordinates.txt');
coordinates = textscan(fd,'%f','delimiter',';,:');
fclose(fd);

% reading product locations
fd = fopen('mazes/tsp products.txt');
itemNumber = textscan(fd,'%f','delimiter',';:,');
fclose(fd);

locations(1,:) = [coordinates{1}(1) coordinates{1}(2)];

for i=3:3:length(itemNumber{1,1})-1
    locations(i/3+1,1)=(itemNumber{1,1}(i));
end

for i=4:3:length(itemNumber{1,1})
    locations((i-1)/3+1,2)=((itemNumber{1,1}(i)));
end

locations(length(locations)+1,:) = [coordinates{1}(3) coordinates{1}(4)];

nodes = length(locations(:,1));

if exist('Pathwayf','var')
    for i = 1:length(Pathwayf)
        for j = 1:length(Pathwayf)
            connectionsf(i,j)=length(Pathwayf{i,j});
        end
    end
end


%% Calculating pathlengths between two points in a maze

% initializing

connection = (nodes*(nodes-1))/2;
if ~exist('connectionsf','var')
  connectionsf = ones(nodes,nodes)*10000;
end
f=@shortestPath;
counter = 0;

% Calculating paths and placing in connections
for x=1:connection
    tempnodes = nodes;
    while(tempnodes > x)
        
            [connections(x,tempnodes),Pathway{x,tempnodes}] = f(locations(x,:),locations(tempnodes,:),runs,iters);

        if tempnodes~=18
             
            if connections(x,tempnodes)<connectionsf(x,tempnodes) || connectionsf(x,tempnodes) == 0
                
                connectionsf(x,tempnodes) = connections(x,tempnodes);
                Pathwayf{x,tempnodes} = Pathway{x,tempnodes};
                connectionsf(tempnodes,x) = connectionsf(x,tempnodes);
                Pathwayf{tempnodes,x} = fliplr(Pathwayf{x,tempnodes});
            end
        end
        
        if ~exist('Pathwayf','var') || tempnodes == 18
            
            [connectionsf(x,tempnodes),Pathwayf{x,tempnodes}] = f(locations(x,:),locations(tempnodes,:),runs,iters);
            connectionsf(tempnodes,x) = connectionsf(x,tempnodes);
            Pathwayf{tempnodes,x} = fliplr(Pathwayf{x,tempnodes});
            
        end
        
        disp(['From: ' num2str(x) ' to '  num2str(tempnodes) ' | PathLength: ' num2str(connectionsf(x,tempnodes))])

        tempnodes = tempnodes-1;
        counter = counter + 1;
    end
    
    
    
end

save('Pathwayf','Pathwayf');

bestPath = 100000000;
for reiter = 1:Reiterate
    %% Creating random chromosomes
    amountchromo = 50;
    chromo = zeros(amountchromo,nodes,1);
    random = 1;
    s = 0;

    %
    for i=1:length(chromo(:,1,1))
        for j=1:length(chromo(1,:,1))
            
            if j>1
                random = round(rand()*18+0.5);
                while sum(random == chromo(i,:,1)) ~= 0
                    random = round(rand()*18+0.5);
                end
                chromo(i,j,1) = random;
            else
                chromo(i,j,1) = 1;
            end
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
                fitness(i,iter) = fitness(i,iter) + connectionsf(point(2),point(1));
                Path{i,j,iter} = Pathwayf{point(1),point(2)};

            end
        end
        
        % first point should be last point
        for i = 1:length(chromo(:,1,iter))
            Path{i,length(chromo(1,:,iter)),iter}=Pathwayf{(chromo(i,end,iter)),chromo(i,1,iter)};
            fitness(i,iter) = fitness(i,iter) + connectionsf((chromo(i,end,iter)),chromo(i,1,iter));
        end

        [pathlength,chromosome]=min(fitness(:,iter));
        fitness(:,iter) = 1./fitness(:,iter)*100000;
        fitnessRatio(:,iter) = fitness(:,iter)./sum(fitness(:,iter));
        
        distribution(1,iter)=0;
        distribution(2:amountchromo+1,iter) = cumsum(fitnessRatio(:,iter));

        amountselections = 2;

        for i=1:amountchromo
           %% selecting
           % selection the best solution - elitism 
           if i==1
                [elit, loca] = max(fitnessRatio(:,iter));
                chromo(i,:,iter+1)=chromo(loca,:,iter);

           end    
               
           if i>1
                for k=1:amountchromo;
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

                    chromo(i,:,iter+1)=tempvect;

                end
    %             if random2 < pcrossover
    % 
    %                 tempvect=zeros(1,nodes);
    % 
    %                 tempvect(1:random)=chromo(select(2),1:random,iter);
    % 
    %                 newchromoloop=random;
    %                 for loopchromo=1:nodes
    % 
    %                     if sum(chromo(select(1),loopchromo,iter) ==  tempvect) == 0;
    %                         newchromoloop=newchromoloop+1;
    %                         tempvect(newchromoloop)= chromo(select(1),loopchromo,iter);
    % 
    %                     end
    %                 end
    % 
    %                 chromo(i+1,:,iter+1)=tempvect;
    % 
    %             end
    % 

                if random2 >= pcrossover
                    chromo(i,:,iter+1)=chromo(select(1),:,iter);
                end
            end
        end




        %% mutation
        pmutation = 0.01;
        for i=1:amountchromo
            for j=1:nodes
                
                if j>1 && j<nodes
                    random = rand;
                else
                    random =2;
                end
                    if random < pmutation
                        change = chromo(i,j,iter+1);
                        random2=round(rand()*18+0.5);
                        chromo(i,j,iter+1)=chromo(i,random2,iter+1);
                        chromo(i,random2,iter+1)=change;

                    end
                
                    
            end
        end
        
        %% andere dingen
        drawnow
        disp(['pathlength: ' num2str(bestPath)]);
        if pathlength < bestPath

            bestPath = pathlength;
            %error(['pathlength: ' num2str(bestPath)]);

            BestPathWay = Path(chromosome,:,iter);

        else
            s=s+1;
            if s>500
                break
            end
            
        end
        
    end
    toc
end
%%
lengte = 0;
for i=1:length(BestPathWay)
    lengte = length(BestPathWay{i})+lengte-1;
end

winddir{1}=num2str(lengte);
winddir{2}=[num2str(BestPathWay{1}(2,1)-1) ', ' num2str(BestPathWay{1}(1,1)-1) ';'];

k=2;
for j = 1:nodes
    for i = 1:length(BestPathWay{j})-1
          
        winddirr = BestPathWay{j}(:,i+1)'-BestPathWay{j}(:,i)';
       
        if i == 1
                        
            k=k+1;
            for p = 1:length(locations)
                if sum(fliplr(locations(p,:))'+1==BestPathWay{j}(:,i))==2
                    stringproduct = ['takeproduct #' num2str(p) ';'];
                    winddir{k} = stringproduct;

                end
            end
        end
        
        
        k=k+1;
        if winddirr(1,1) == 0 && winddirr(1,2) == -1
            winddir{k} = '2;';
        
        
        elseif winddirr(1,1) == 0 && winddirr(1,2) == 1
            winddir{k} = '0;';
        
        
        elseif winddirr(1,1) == 1 && winddirr(1,2) == 0
             winddir{k} = '3;';
        
        
        elseif winddirr(1,1) == -1 && winddirr(1,2) == 0
             winddir{k} = '1;';
        
        end
        
    end
end

fid=fopen('5_actions_TSP.txt','w');


for i=1:length(winddir);
    fprintf(fid,'%s\n',winddir{i});
end

fclose(fid);


