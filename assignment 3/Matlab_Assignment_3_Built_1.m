%% initializing
clear all
close all 
clc

%% Input Maze

% input of the mazes 
easyMaze = dlmread('Mazes\easy maze.txt','',1,0);
mediumMaze = dlmread('Mazes\medium maze.txt','',1,0);
hardMaze = dlmread('Mazes\hard maze.txt','',1,0);
training1Maze = dlmread('Mazes\training\training1.txt','',1,0);

% input of the starting and ending coordinates of the maze
easyMaze_c = load('Mazes\easy coordinates.txt');
mediumMaze_c = load('Mazes\medium coordinates.txt');
hardMaze_c = load('Mazes\hard coordinates.txt');
training1Maze_c = load('Mazes\Training\training1 coordinates.txt');

%% Variables

% input
iterations = 1;
ants = 1;
pheromone = 100;
evaporation = 0.1;
stopcriterian = 0.90;
windDir(1,:)=[0,1];
windDir(2,:)=[-1,0];
windDir(3,:)=[0,-1];
windDir(4,:)=[1,0];

% calculation
beginPos = training1Maze_c(1,:)+1;
endPos = training1Maze_c(2,:)+1;

mazeSize = size(training1Maze);

pherMatrix=training1Maze;

%% For now: 1 ant walking randomly through maze

Pos = beginPos;
q=1;

% While ants are running towards end
while Pos(1,1) ~= endPos(1,1) || Pos(1,2) ~= endPos(1,2)
    
    % Number of iterations is higher
    q=q+1;
    
    % List of previous position (important for not going to previous place)
    PosList(1,:) = beginPos;
    
    
    % For every possible next position
    for i=1:4
        
        
        % Calculate the position next the the current position
        choices(i,:)= Pos+windDir(i,:);
        
        
        % If is it a valid position
        if  choices(i,1) >= 1 && ...
            choices(i,1) <= mazeSize(1,1) && ...
            choices(i,2) >= 1 && ...
            choices(i,2) <= mazeSize(1,2)
            
            % insert the property (wall or road)(0 or 1)
            choicesProperty(i)=training1Maze(choices(i,1),choices(i,2));

        else
            
            % if it is not valid, it should be 0
            choicesProperty(i)=0;

        end
        
    end
    
    
    %for all choices
    for i=1:length(choicesProperty)
        
        
        % If the choice is 1
        if choicesProperty(i) == 1;
 
            %calculate the newPos and check if it is already visited
            posNew=Pos+windDir(i,:);
  
            I = sum(PosList(:, 1) == posNew(1) & PosList(:, 2) == posNew(2));
            
            % If so, then change its choiceProperty to 0
            if I 
                
                choicesProperty(i)=2;
                
            end
                
        end
        
        
        
    end
    
    % creating a matrix for the pheromones in all posible directions
    pher = zeros(1,4);
    
    % loop through the directions
    for i=1:length(choicesProperty)        
        
        % if a direction is a posibility 
        if choicesProperty(i) == 1;
            
            % calculate the coordinates of the possible point
            NewPos=Pos + windDir(i,:);
            
            % determine how much pheremones are on this spot
            pher(:,i) = pherMatrix(NewPos(1,1),NewPos(1,2));
            
        end
    end
    
    % creating the probability's of choosing the directions
    probs = pher/sum(pher);
    
    % creating a random number
    random = rand();
    
    % if the random number lies between 0 and the probability of first direcition
    % choose that one elseif the random number lies between the probability
    % of the first and second direction choose that one etc
    
    if random<probs(1)
        Pos=Pos+windDir(1,:);
    elseif random>probs(1) && rand<probs(1)+probs(2)
        Pos=Pos+windDir(2,:);
    elseif random>probs(1)+probs(2) && rand<probs(1)+probs(2)+probs(3)
        Pos=Pos+windDir(3,:);
    elseif random>probs(1)+probs(2)+probs(3)
        Pos=Pos+windDir(4,:);
    end

   % keep track of the route the ant took
   PosList(q,:) = Pos;
    
end

%% plotting the route of an ant and plotting the maze with begin position as well as end position

% plot the route of the ant
plot(PosList(:,2),-PosList(:,1)+1,'r')

% setting the axes ranges
axis([-1 mazeSize(1,2)+2 -mazeSize(1,1)-1 2])
hold on

% looping through the whole maze
for i=0:mazeSize(1,1)+1;
    for j=0:mazeSize(1,2)+1;
        
        % if the i or j is not in the matrix create an X
        if i < 1 || i > mazeSize(1,1) || j < 1 || j > mazeSize(1,2)
            
            plot(j,-i+1,'x','linewidth',20)

        
        else
            
            % if the position on the maze has a wall (0) plot an X there
            if training1Maze(i,j)==0
                
                plot(j,-i+1,'x','linewidth',20)
                
            % if  the postion on the maze is the end position plot a green
            % square
            elseif i== endPos(1,1) && j== endPos(1,2)
                
                plot(j,-i+1, 'gs','linewidth',5)
            
            % if the position on the maze is the begin position plot a
            % magenta square
            elseif i== beginPos(1,1) && j== beginPos(1,2)
                
                plot(j,-i+1, 'ms','linewidth',5)
                
            end
            
        end
    end
    
    % making plot area square
    axis square
end