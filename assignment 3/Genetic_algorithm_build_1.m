clear all;
fd = fopen('mazes/tsp products.txt');
itemNumber = textscan(fd,'%f','delimiter',';:,');
fclose(fd);

for i=3:3:length(itemNumber{1,1})-1
    locations(i/3,1)=(itemNumber{1,1}(i));
end

for i=4:3:length(itemNumber{1,1})
    locations((i-1)/3,2)=((itemNumber{1,1}(i)));
end

%% Calculating pathlengths between two points in a maze

% initializing
nodes = length(locations(:,1));
connection = (nodes*(nodes-1))/2;
connections = zeros(1,nodes);
f=@shortestPath;
counter = 0;

% Calculating paths and placing in connections
for x=1:connection
    tempnodes = nodes;
    while(tempnodes > x)
        connections(x,tempnodes) = f(locations(x,:),locations(tempnodes,:),1,6);
        connections(tempnodes,x) = connections(x,tempnodes);

        tempnodes = tempnodes-1;
        counter = counter + 1;
        disp(counter)
    end
end

%% Creating random chromosomes
amountchromo = 10;
chromo = zeros(amountchromo,nodes);
random = 1;

%
for i=1:length(chromo(:,1))
    for j=1:length(chromo(1,:))
        random = round(rand()*18+0.5);
        while sum(random == chromo(i,:)) ~= 0
            random = round(rand()*18+0.5);
        end
        chromo(i,j) = random;
    end
end


%% fitness ratio
fitness = zeros(amountchromo,1);

for i=1:length(chromo(:,1))
    for j=1:length(chromo(1,:))-1
        point(1) = chromo(i,j);
        point(2) = chromo(i,j+1);
        fitness(i) = fitness(i) + connections(point(2),point(1));
    end
end
fitness = 1./fitness;
fitnessRatio = fitness./sum(fitness);

%% selecting
% choosing the new population
distribution(1)=0;
distribution(2:amountchromo+1) = cumsum(fitnessRatio);
amountselections = 2;
select =zeros(1,amountselections);

while length(unique(select)) ~= amountselections
    for i=1:amountselections
        random = rand;
        for j=1:amountchromo
            if distribution(j) < random && distribution(j+1) > random
                select(i) = j;
                break
            end
        end
    end
end

%% mutation
pmutation = 0.005
for i=1:length(select)
    for j=1:nodes
        random = rand;
        if rand < pmutation
end


