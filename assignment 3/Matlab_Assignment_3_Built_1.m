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

pherMarix=training1Maze;

%%

Pos = beginPos;
q=0;

% While ants are running towards end
while Pos ~= endPos
    
    % Number of iterations is higher
    q=q+1;
    
    % List of previous position (important for not going to previous place)
    PosList(q,:) = Pos;
    
    
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
    
    for i=1:length(choicesProperty)
        counter = 0;
        pher = zeros(1,4);
        if choiceProperty(i) == 1;
            pheromone(i) = pherMatrix(Pos + windDir(i,:));
            counter = counter+1;
        end
    end
    
    
        
    
    
   %TIJDELIJK !!!! 
   Pos=endPos;
    
end
