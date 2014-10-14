%% initializing
clear all
close all 
clc

%% Input Maze

% input of the mazes 
Maze = dlmread('Mazes\medium maze.txt','',1,0);

% input of the starting and ending coordinates of the maze
Maze_c = load('Mazes\medium coordinates.txt');

%% Variables

% input
iterations = 1000;
ants = 10;
pheromone = 30;
evaporation = 0.95;
stopcriterian = 0.90;
windDir(1,:)=[0,1];
windDir(2,:)=[-1,0];
windDir(3,:)=[0,-1];
windDir(4,:)=[1,0];

% calculation
beginPos = Maze_c(1,:)+1;
endPos = Maze_c(2,:)+1;

mazeSize = size(Maze);

pherMatrix=Maze;

%% For now: 1 ant walking randomly through maze

% Creating matrices
Pos = cell(iterations,ants);
PosList = cell(iterations,ants);
PosList2 = cell(iterations,ants);
choices = zeros(4,2);
choicesProperty = zeros(1,4);
pherCurMatrix=cell(iterations,ants);
posfinalList = zeros(2,100000);

for itt = 1:iterations
    
%     clf()
%     ants = round(ants/sqrt(itt)) + 50;
%     pheromone = itt *pheromone;
    
    for ant = 1:ants
        r=1;
        Pos{itt,ant} = beginPos;
        q=1;
        PosList{itt,ant}(:,1) = beginPos;
        PosList2{itt,ant}=[];
        
        % While ants are running towards end
        while Pos{itt,ant}(1,2) ~= endPos(1,1) || Pos{itt,ant}(1,1) ~= endPos(1,2)

            % Number of iterations is higher
            q=q+1;
            
            % For every possible next position
            for i=1:4
                
                % Calculate the position next the the current position
                choices(i,:)= Pos{itt,ant}+windDir(i,:);
                
                % If is it a valid position
                if  choices(i,1) >= 1 && ...
                    choices(i,1) <= mazeSize(1,1) && ...
                    choices(i,2) >= 1 && ...
                    choices(i,2) <= mazeSize(1,2)

                    % insert the property (wall or road)(0 or 1)
                    choicesProperty(i)=Maze(choices(i,1),choices(i,2));
                    
                        % If the choice is 1
                    if choicesProperty(i) == 1;

                        %calculate the newPos and check if it is already visited
                        posNew=Pos{itt,ant}+windDir(i,:);

                        % If so, then change its choiceProperty to
                        if sum(PosList{itt,ant}(1, :) == posNew(1) & PosList{itt,ant}(2, :) == posNew(2));

                            choicesProperty(i)=2;

                        end

                    end
                    
                else

                    % if it is not valid, it should be 0
                    choicesProperty(i)=0;

                end

            end



                % In case the number of possible unchecked routes is zero

            p=1;
            while sum(choicesProperty == 1) == 0

                Pos{itt,ant}=PosList{itt,ant}(:,q-p)';
                PosList2{itt,ant}(:,r) = Pos{itt,ant};
                    % For every possible next position
                for i=1:4


                    % Calculate the position next the the current position
                    choices(i,:)= Pos{itt,ant}+windDir(i,:);


                    % If is it a valid position
                    if  choices(i,1) >= 1 && ...
                        choices(i,1) <= mazeSize(1,1) && ...
                        choices(i,2) >= 1 && ...
                        choices(i,2) <= mazeSize(1,2)

                        % insert the property (wall or road)(0 or 1)
                        choicesProperty(i)=Maze(choices(i,1),choices(i,2));
                    
                         % If the choice is 1
                    if choicesProperty(i) == 1;

                        %calculate the newPos and check if it is already visited
                        posNew=Pos{itt,ant}+windDir(i,:);

                        % If so, then change its choiceProperty to 2
                        if sum(PosList{itt,ant}(1, :) == posNew(1) & PosList{itt,ant}(2, :) == posNew(2));

                            choicesProperty(i)=2;

                        end

                    end
                        
                    else

                        % if it is not valid, it should be 0
                        choicesProperty(i)=0;

                    end  

                end

                r=r+1;
                p=p+1;

            end
            
                        % creating a matrix for the pheromones in all posible directions
            pher = zeros(1,4);

            % In case the number of possible unchecked routes is at least 1

            if sum(choicesProperty == 1) > 0
                % loop through the directions
                for i=1:length(choicesProperty)        

                    % if a direction is a posibility 
                    if choicesProperty(i) == 1;

                        % calculate the coordinates of the possible point
                        NewPos=Pos{itt,ant} + windDir(i,:);

                        % determine how much pheremones are on this spot
                        pher(:,i) = pherMatrix(NewPos(1,1),NewPos(1,2));

                    end

                end
            end


            % creating the probability's of choosing the directions
            probs{itt,ant}(:,q) = pher./sum(pher);

            % creating a random number
            random = rand();

            % if the random number lies between 0 and the probability of first direcition
            % choose that one elseif the random number lies between the probability
            % of the first and second direction choose that one etc

            if random<=probs{itt,ant}(1,q)
                Pos{itt,ant}=Pos{itt,ant}+windDir(1,:);
            elseif random>probs{itt,ant}(1,q) && random<=probs{itt,ant}(1,q)+probs{itt,ant}(2,q)
                Pos{itt,ant}=Pos{itt,ant}+windDir(2,:);
            elseif random>probs{itt,ant}(1,q)+probs{itt,ant}(2,q) && random<=probs{itt,ant}(1,q)+probs{itt,ant}(2,q)+probs{itt,ant}(3,q)
                Pos{itt,ant}=Pos{itt,ant}+windDir(3,:);
            elseif random>probs{itt,ant}(1,q)+probs{itt,ant}(2,q)+probs{itt,ant}(3,q)
                Pos{itt,ant}=Pos{itt,ant}+windDir(4,:);
            end

           % keep track of the route the ant took
           PosList{itt,ant}(:,q) = Pos{itt,ant};


           if length(PosList{itt,ant})>100000
               break
           end

        end
        
        pher_pp = pheromone/(length(PosList{itt,ant}));
        
        pherCurMatrix{itt,ant} = zeros(mazeSize(1,1),mazeSize(1,2));
        
        for i=1:length(PosList{itt,ant})
            
            pherCurMatrix{itt,ant}(PosList{itt,ant}(1,i),PosList{itt,ant}(2,i)) = pherCurMatrix{itt,ant}(PosList{itt,ant}(1,i),PosList{itt,ant}(2,i))+pher_pp;

        end 
        
        PosList{itt,ant}(:,length(PosList{itt,ant})+1:length(PosList2{itt,ant})+length(PosList{itt,ant})) = PosList2{itt,ant};    

        if length(PosList{itt,ant}) < length(posfinalList)
            
            posfinalList = PosList{itt,ant};
            a = ant;
            b = itt;
            
        end
            
%         disp(ant)

    end
    

    
    for ant = 1:ants
        
        pherMatrix = pherMatrix + pherCurMatrix{itt,ant};
            
    end
    
    c = pherCurMatrix{itt,ant} == 0;
    c=c*evaporation;
    d = c == 0;
    c=c+d;
    
    pherMatrix = pherMatrix.*c;
    
    disp(itt)
    
    maxpher = max(max(pherMatrix));
minpher = min(min(pherMatrix));

% for i = 1:mazeSize(1,1)
%     for j = 1:mazeSize(1,2)
%         
%         plot(j,-i,'x','Color',[1-pherMatrix(i,j)/maxpher,((pherMatrix(i,j)/maxpher)),0],'LineWidth',5)
%         hold on
%         
%     end
% end
% 
% drawnow
end
%% plotting theroute of an ant and plotting the maze with begin position as well as end position
    
    % plot the route of the ant
    plot(PosList(2,:),-PosList(1,:)+1,'Color',[rand(),rand(),rand()])
    hold on

    % setting the axes ranges
    axis([-1 mazeSize(1,2)+2 -mazeSize(1,1)-1 2])


% looping through the whole maze
for i=0:mazeSize(1,1)+1;
    for j=0:mazeSize(1,2)+1;
        
        % if the i or j is not in the matrix create an X
        if i < 1 || i > mazeSize(1,1) || j < 1 || j > mazeSize(1,2)
            
            plot(j,-i+1,'x','linewidth',mazeSize(1,1)/10)

        
        else
            
            % if the position on the maze has a wall (0) plot an X there
            if Maze(i,j)==0
                
                plot(j,-i+1,'x','linewidth',mazeSize(1,2)/10)
                
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
    

end

    hFig = figure(1);
    set(hFig, 'Position', [0 0 mazeSize(1,2)*15 mazeSize(1,1)*15]);
    movegui(hFig,'center')
    % making plot area square
%     axis square

close(figure(1))
%% Debug
 h = figure(2);
 hold on
 hFig = figure(2);
 set(hFig, 'Position', [0 0 mazeSize(1,2)*15 mazeSize(1,1)*15]);
 movegui(hFig,'center')
 
 % setting the axes ranges



    
   
    % setting the axes ranges
    axis([-1 mazeSize(1,2)+2 -mazeSize(1,1)-1 2])

    % looping through the whole maze
    for i=0:mazeSize(1,1)+1;
        for j=0:mazeSize(1,2)+1;

            % if the i or j is not in the matrix create an X
            if i < 1 || i > mazeSize(1,1) || j < 1 || j > mazeSize(1,2)

                plot(j,-i+1,'x','linewidth',mazeSize(1,1)/10)
                hold on
                axis([-1 mazeSize(1,2)+2 -mazeSize(1,1)-1 2])
                
            else

                % if the position on the maze has a wall (0) plot an X there
                if Maze(i,j)==0

                    plot(j,-i+1,'x','linewidth',mazeSize(1,2)/10)

                end

            end
        end

    end
    
Posfinal=beginPos;
q=1;
posfinalList(:,q) = Posfinal;

choicesfinal = zeros(1,4);
for i = 1 : length(posfinalList)-length(PosList2{b,a})
    
    % plot the route of the ant
    ant = plot(posfinalList(2,i),-posfinalList(1,i)+1,'xr','LineWidth',5);
    plot(posfinalList(2,i),-posfinalList(1,i)+1,'xm','LineWidth',2)
    hold on
    
    drawnow
    delete(ant)
    
end


%%


