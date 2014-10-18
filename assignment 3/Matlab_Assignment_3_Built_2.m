%% initializing
clear all
close all 
clc

%% Input Maze

% input of the mazes 
Maze = dlmread('Mazes\hard maze.txt','',1,0);

% input of the starting and ending coordinates of the maze
Maze_c = load('Mazes\hard coordinates.txt');

%% Variables

% input
iterations = 100;
pheromone = 850;
evaporation = 0.5;
windDir(1,:)=[0,1];
windDir(2,:)=[-1,0];
windDir(3,:)=[0,-1];
windDir(4,:)=[1,0];

posfinalList = zeros(2,1000000);

k=0;
beginPos = Maze_c(1,:)+1;
endPos = Maze_c(2,:)+1;
while 1==1
    mazeSize = size(Maze);
    k=k+1;
    PosBestOwn = zeros(2,1000000);
    s=0;
    pherMatrix=Maze;
    ants=4;

    clear PosList
    clear PosList2
    clear pherCurMatrix
    clear LastDifChoice

    %% ants movining through maze
    % 
    % Creating matrices
    PosList = cell(1,1);
    PosList2 = cell(1,1);
    choices = zeros(4,2);
    pherCurMatrix=cell(iterations);

    for i = 1:iterations
        for j = 1:ants
            PosList{i,j} = zeros(2,1);
            PosList2{i,j} = zeros(2,1);
            pherCurMatrix{i,j} = zeros(mazeSize(1,1),mazeSize(1,2));

            drawnow
        end

    end
    itt=0;
    o=0;
    while 1==1
        itt=itt+1;
        o=o+1;
        if o>40 && ants > 1
            ants=ants-1;
            o=0;
        end
    %     ants = round(ants/sqrt(itt)) + 50;
    %     pheromone = itt *pheromone;
    %     
        for ant = 1:ants

            q=1;
            r=1;

            Pos = beginPos;
            PosList{itt,ant}(:,1) = beginPos;


            % While ants are running towards end
            while Pos(1,2) ~= endPos(1,1) || Pos(1,1) ~= endPos(1,2)
                tic
                % Number of iterations is higher
                q=q+1;


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
                        choicesProperty(i,q)=Maze(choices(i,1),choices(i,2));

                            % If the choice is 1
                        if choicesProperty(i,q) == 1;

                           % If so, then change its choiceProperty to 2
                            if sum(PosList{itt,ant}(1, :) == choices(i,1) & PosList{itt,ant}(2, :) == choices(i,2));

                                choicesProperty(i,q)=2;

                            end

                        end

                    else

                        % if it is not valid, it should be 0
                        choicesProperty(i,q)=0;

                    end

                end



                % In case the number of possible unchecked routes is zero

                p=1;
                while sum(choicesProperty(:,q) == 1) == 0

                    Pos=PosList{itt,ant}(:,q-p)';
                    PosList2{itt,ant}(:,r) = Pos;
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
                            choicesProperty(i,q)=Maze(choices(i,1),choices(i,2));

                             % If the choice is 1
                        if choicesProperty(i,q) == 1;


                            % If so, then change its choiceProperty to 2
                            if sum(PosList{itt,ant}(1, :) == choices(i,1) & PosList{itt,ant}(2, :) == choices(i,2));

                                choicesProperty(i,q)=2;

                            end

                        end


                        else

                            % if it is not valid, it should be 0
                            choicesProperty(i,q)=0;

                        end  

                    end

                    r=r+1;
                    p=p+1;

                end

                            % creating a matrix for the pheromones in all posible directions
                pher = zeros(1,4);



                % In case the number of possible unchecked routes is at least 1

                if sum(choicesProperty(:,q) == 1) > 0
                    % loop through the directions
                    for i=1:length(choicesProperty(:,q))        

                        % if a direction is a posibility 
                        if choicesProperty(i,q) == 1;

                            % calculate the coordinates of the possible point
                            NewPos=Pos + windDir(i,:);

                            % determine how much pheremones are on this spot
                            pher(:,i) = pherMatrix(NewPos(1,1),NewPos(1,2));

                        end

                    end
                end
                
               

                if sum((choicesProperty(:,q-1) == 1) + (choicesProperty(:,q-1) == 2))>=3 && q>2 
                    
                    if choice(q-1)~=choice(q-2)
                       LastDifChoice = choice(q-2);
                    end
                    
                    if choice(q-1)==choice(q-2) && sum((choicesProperty(:,q-1) == 1) + (choicesProperty(:,q-1) == 2))==3 && sum((choicesProperty(:,q) == 1) + (choicesProperty(:,q) == 2))==4
                        
                        pher = pher .* (choicesProperty(:,q-1) == 0)';
                    
                    end
                    
                    if exist('LastDifChoice','var') 
                        if (LastDifChoice == 1 || LastDifChoice == 2) 
                            pher(LastDifChoice+2) = 0;
                        elseif LastDifChoice == 3
                            pher(1) = 0;
                        elseif LastDifChoice == 4
                            pher(2) = 0;
                        end
                    end
                    

               
                        
                end

                if sum(pher)~=0

                    % creating the probability's of choosing the directions
                    probs{itt,ant}(:,q) = pher./sum(pher);

                else

                    probs{itt,ant}(:,q) = choicesProperty(:,q) == 1;
                end


                % creating a random number
                random = rand();


                % if the random number lies between 0 and the probability of first direcition
                % choose that one elseif the random number lies between the probability
                % of the first and second direction choose that one etc

                if random<=probs{itt,ant}(1,q)
                    Pos=Pos+windDir(1,:);
                    choice(q)=1;
                elseif random>probs{itt,ant}(1,q) && random<=probs{itt,ant}(1,q)+probs{itt,ant}(2,q)
                    Pos=Pos+windDir(2,:);
                    choice(q)=2;
                elseif random>probs{itt,ant}(1,q)+probs{itt,ant}(2,q) && random<=probs{itt,ant}(1,q)+probs{itt,ant}(2,q)+probs{itt,ant}(3,q)
                    Pos=Pos+windDir(3,:);
                    choice(q)=3;
                elseif random>probs{itt,ant}(1,q)+probs{itt,ant}(2,q)+probs{itt,ant}(3,q)
                    Pos=Pos+windDir(4,:);
                    choice(q)=4;
                end


               % keep track of the route the ant took
               PosList{itt,ant}(:,q) = Pos;


               if length(PosList{itt,ant})>10000
                   break
               end

            end


            pher_pp = pheromone/(length(PosList{itt,ant}));

            pherCurMatrix{itt,ant} = zeros(mazeSize(1,1),mazeSize(1,2));

            for i=1:length(PosList{itt,ant})

                pherCurMatrix{itt,ant}(PosList{itt,ant}(1,i),PosList{itt,ant}(2,i)) = pherCurMatrix{itt,ant}(PosList{itt,ant}(1,i),PosList{itt,ant}(2,i))+pher_pp;

            end
            PosList2{itt,ant} = [] ;
            PosList{itt,ant}(:,length(PosList{itt,ant})+1:length(PosList2{itt,ant})+length(PosList{itt,ant})) = PosList2{itt,ant};    
            s=s+1;
            if length(PosList{itt,ant}) < length(posfinalList)

                posfinalList = PosList{itt,ant};
                a = ant;
                b = itt;
                s=0;

            end

            if length(PosList{itt,ant})< length(PosBestOwn)
                PosBestOwn = PosList{itt,ant};
                s=0;
            end



        disp(['iteration: ' num2str(itt) ' | ant: ' num2str(ant) ' | # no better solution: ' num2str(s) ' | Best pathlength: ' num2str(length(posfinalList)) ' | Current best: ' num2str(length(PosBestOwn)) ' | Current Pathlength: ' num2str(length(PosList{itt,ant}))]);
        drawnow
        end

        for antt = 1:ants

            pherMatrix = (pherMatrix + pherCurMatrix{itt,ant}).*(pherMatrix>10^-20);

        end

        c = pherCurMatrix{itt,ant} == 0;
        c=c*(1-evaporation);
        d = c == 0;
        c=c+d;

        pherMatrix = pherMatrix.*c;

        drawnow

        maxpher = max(max(pherMatrix));
        minpher = min(min(pherMatrix));

        if itt>3
            if s>ants*1000 || (length(PosList{itt,1}) == length(PosList{itt-1,1}) && length(PosList{itt,1}) == length(PosList{itt-2,1}) && length(PosList{itt,1}) == length(PosList{itt-3,1}));
                converged(k) = length(PosList{itt,ant});
                break
            end
        end
    end
end
%% plotting theroute of an ant and plotting the maze with begin position as well as end position
    
%     % plot the route of the ant
%     plot(PosList(2,:),-PosList(1,:)+1,'Color',[rand(),rand(),rand()])
%     hold on
% 
%     % setting the axes ranges
%     axis([-1 mazeSize(1,2)+2 -mazeSize(1,1)-1 2])


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

 set(h, 'Position', [0 0 mazeSize(1,2)*15 mazeSize(1,1)*15]);
 movegui(h,'center')
 
% setting the axes ranges
axis([-1 mazeSize(1,2)+2 -mazeSize(1,1)-1 2])

    % looping through the whole maze
    f=0;
    MatrixPlot=zeros(2,mazeSize(1,1)*mazeSize(1,2));
    for i=0:mazeSize(1,1)+1;
        for j=0:mazeSize(1,2)+1;
            
            % if the i or j is not in the matrix create an X
            if i < 1 || i > mazeSize(1,1) || j < 1 || j > mazeSize(1,2)
                f=f+1;
                MatrixPlot(:,f) = [j ; -i+1];
                
            else

                % if the position on the maze has a wall (0) plot an X there
                if Maze(i,j)==0
                    f=f+1;
                    MatrixPlot(:,f) = [j ; -i+1];

                end

            end
        end
       
    end
    
    h=plot(MatrixPlot(1,1:1),MatrixPlot(2,1:1),'x','linewidth',mazeSize(1,2)/10);
    for f=2:mazeSize(1,2):f
        set(h,'Xdata',MatrixPlot(1,1:f+mazeSize(1,2)),'Ydata',MatrixPlot(2,1:f+mazeSize(1:2)));
        drawnow
    end

Posfinal=beginPos;
q=1;
posfinalList(:,q) = Posfinal;

choicesfinal = zeros(1,4);
i=1;
h=plot(posfinalList(2,i),-posfinalList(1,i)+1,'xr','LineWidth',5);
j=plot(posfinalList(2,i),-posfinalList(1,i)+1,'xm','LineWidth',2);

while ~sum(posfinalList(1, i) == endPos(1,2) & posfinalList(2, i) == endPos(1,1));
    % plot the route of the ant
    
    set(h,'Xdata',posfinalList(2,i),'Ydata',-posfinalList(1,i)+1);
    set(j,'Xdata',posfinalList(2,1:i),'Ydata',-posfinalList(1,1:i)+1);

    hold on
    
    drawnow
    i=i+1;

end


%%

 h = figure(2);
 hold on
 set(h, 'Position', [0 0 mazeSize(1,2)*15 mazeSize(1,1)*15]);
 movegui(h,'center')

for i = 1:mazeSize(1,1)
    
    for j = 1:mazeSize(1,2)
        
        plot(j,-i,'x','Color',[1-pherMatrix(i,j)/maxpher,((pherMatrix(i,j)/maxpher)),0],'LineWidth',18)
        hold on
        
    end
end

drawnow