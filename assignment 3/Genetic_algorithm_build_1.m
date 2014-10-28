tic
clear all;
clc

%% loading previous point to point route data

% load the previous point to point route data
fid = fopen('Pathwayf.mat');
if fid~=-1
    load('Pathwayf.mat','Pathwayf');
end

% if it existed create a matrix with the point to point route lengths in it
if exist('Pathwayf','var')
    for i = 1:length(Pathwayf)
        for j = 1:length(Pathwayf)
            connectionsf(i,j)=length(Pathwayf{i,j});
        end
    end
end

% if it did not exist create a new matrix for the point to point route lengths
if ~exist('connectionsf','var')
  connectionsf = ones(nodes,nodes)*10000;
end

%% list of variables
iterations = 1000;              % number of iterations for the tsp
Reiterate = 20;                  % number of runs for the tsp
runs = 100;                     % number of runs for the point to point routes
iters = 40;                     % number of iterations for the point to point routes
f = @shortestPath;              % function for calculating the fastest route between to points

%% reading in starting, end and products locations

% start and end coordinates
fd = fopen('mazes/hard coordinates.txt');
coordinates = textscan(fd,'%f','delimiter',';,:');
fclose(fd);

% reading product locations
fd = fopen('mazes/tsp products.txt');
itemNumber = textscan(fd,'%f','delimiter',';:,');
fclose(fd);

% writing starting point to locations
locations(1,:) = [coordinates{1}(1) coordinates{1}(2)];

% writing product x coordinates to locations
for i=3:3:length(itemNumber{1,1})-1
    locations(i/3+1,1)=(itemNumber{1,1}(i));
end

% writing product y coordinates to locations
for i=4:3:length(itemNumber{1,1})
    locations((i-1)/3+1,2)=((itemNumber{1,1}(i)));
end

% writing end point to locations
locations(length(locations)+1,:) = [coordinates{1}(3) coordinates{1}(4)];

% calculating the ammount of nodes to be taken into account
nodes = length(locations(:,1));

% calculating number of routes possible bewteen points
connection = (nodes*(nodes-1))/2;



% %% Calculating pathlengths between two points in a maze
% 
% % looping through all possible routes between to points
% for x=1:connection
%     
%     tempnodes = nodes;
%     while(tempnodes > x)
%         
%             [connections(x,tempnodes),Pathway{x,tempnodes}] = f(locations(x,:),locations(tempnodes,:),runs,iters);
% 
%         if tempnodes~=nodes
%              
%             if connections(x,tempnodes)<connectionsf(x,tempnodes) || connectionsf(x,tempnodes) == 0
%                 
%                 connectionsf(x,tempnodes) = connections(x,tempnodes);
%                 Pathwayf{x,tempnodes} = Pathway{x,tempnodes};
%                 connectionsf(tempnodes,x) = connectionsf(x,tempnodes);
%                 Pathwayf{tempnodes,x} = fliplr(Pathwayf{x,tempnodes});
%             end
%         end
%         
%         if ~exist('Pathwayf','var') || tempnodes == nodes
%             
%             [connectionsf(x,tempnodes),Pathwayf{x,tempnodes}] = f(locations(x,:),locations(tempnodes,:),runs,iters);
%             connectionsf(tempnodes,x) = connectionsf(x,tempnodes);
%             Pathwayf{tempnodes,x} = fliplr(Pathwayf{x,tempnodes});
%             
%         end
%         
%         disp(['From: ' num2str(x) ' to '  num2str(tempnodes) ' | PathLength: ' num2str(connectionsf(x,tempnodes))])
% 
%         tempnodes = tempnodes-1;
%     end
%     
%     
%     
% end
% 
% save('Pathwayf','Pathwayf');
% 
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
            
            if j>1 && j<nodes
                random = round(rand()*nodes+0.5);
                while sum(random == chromo(i,:,1)) ~= 0 || random == 1 || random == 20
                    random = round(rand()*(nodes)+0.5);
                end
                chromo(i,j,1) = random;
            elseif j == 1
                chromo(i,j,1) = 1;
            elseif j == nodes;
                chromo(i,j,1) = 20;
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
        
%         % first point should be last point
%         for i = 1:length(chromo(:,1,iter))
%             Path{i,length(chromo(1,:,iter)),iter}=Pathwayf{(chromo(i,end,iter)),chromo(i,1,iter)};
%             fitness(i,iter) = fitness(i,iter) + connectionsf((chromo(i,end,iter)),chromo(i,1,iter));
%         end

        [pathlength,chromosome]=min(fitness(:,iter));
        fitness(:,iter) = 1./fitness(:,iter).^2*100000;
        fitnessRatio(:,iter) = fitness(:,iter)./sum(fitness(:,iter));
        
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
                random = round(rand()*nodes+0.5);

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
        pmutation = 0.006;
        for i=1:amountchromo
            for j=1:nodes
                
                if j>1 && j<nodes-1
                    random = rand;
                else
                    random =2;
                end
                
                if random < pmutation
                    
                        change = chromo(i,j,iter+1);
                        random2 = round(rand()*(nodes)+0.5);
                        
                    while random2 == 1 || random2 == 20
                        random2 = round(rand()*(nodes)+0.5);
                    end
                        chromo(i,j,iter+1)=chromo(i,random2,iter+1);
                        chromo(i,random2,iter+1)=change;

                end
                
                    
            end
        end
        
        %% andere dingen
        drawnow
        if pathlength < bestPath

            bestPath = pathlength;
            disp(['pathlength: ' num2str(bestPath)]);

            BestPathWay = Path(chromosome,:,iter);

        else
            s=s+1;
            if s>1000
                break
            end
            
        end
        
    end
    toc
end
%% Writing to file

lengte = 0;
for i=1:length(BestPathWay)
    lengte = length(BestPathWay{i})+lengte-1;
end

lengte = lengte + nodes -2;

winddir{1}=num2str([num2str(lengte) ';']);
winddir{2}=[num2str(BestPathWay{1}(2,1)-1) ', ' num2str(BestPathWay{1}(1,1)-1) ';'];

k=2;
for j = 1:nodes-1
    for i = 1:length(BestPathWay{j})-1
          
        winddirr = BestPathWay{j}(:,i+1)'-BestPathWay{j}(:,i)';

        for p = 1:nodes
            if sum(fliplr(locations(p,:))'+1==BestPathWay{j}(:,i))==2 && p~=1 && p~=nodes
                k=k+1;
                stringproduct = ['take product #' num2str(p-1) ';'];
                winddir{k} = stringproduct;

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
        
        else
            disp(winddirr)
            disp(k)
        end
        
    end
end

fid=fopen('5_actions_TSP.txt','w');


for i=1:length(winddir);
    fprintf(fid,'%s\n',winddir{i});
end

fclose(fid);


