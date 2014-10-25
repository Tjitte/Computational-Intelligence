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

%%
nodes = length(locations(:,1));
connection = (nodes*(nodes-1))/2;
connections = zeros(1,nodes);
f=@shortestPath;
counter = 0;
for x=1:connection
    tempnodes = nodes;
    while(tempnodes > x)
        connections(x,tempnodes) = f(locations(x,:),locations(tempnodes,:),1);
        tempnodes = tempnodes-1;
        counter = counter + 1;
        disp(counter)
    end
end