
%% Input Maze

mazestring = 'Mazes/easy';

% input of the mazes 
Maze = dlmread([mazestring ' maze.txt'],'',1,0);

% input of the starting and ending coordinates of the maze
Maze_c = load([mazestring ' coordinates.txt']);

%% Variables

% input
pheromone = 850;
evaporation = 0.9;
windDir(1,:)=[0,1];
windDir(2,:)=[-1,0];
windDir(3,:)=[0,-1];
windDir(4,:)=[1,0];

% setting initial pathlength to a very high number
pathLength = zeros(2,1000000);

% setting begin position end end position
beginPos = Maze_c(1,:)+1;
endPos = Maze_c(2,:)+1;

k=0;

while 1==1

    Maze = dlmread([mazestring ' maze.txt'],'',1,0);
    backupMaze = dlmread([mazestring ' maze.txt'],'',1,0);
    mazeSize = size(Maze);
    k=k+1;
    PosBestOwn = zeros(2,1000000);
    s=0;
    pherMatrix=Maze;
    ants=2;
    
    
    clear PosList
    clear PosList2
    clear pherCurMatrix
    clear LastDifChoice

    %% ants movining through maze
    % 
    % Creating matrices
    PosList = cell(1,1);
    choices = zeros(4,2);
    pherCurMatrix=cell(100);

    for i = 1:100
        for j = 1:ants
            PosList{i,j} = zeros(2,1);
            pherCurMatrix{i,j} = zeros(mazeSize(1,1),mazeSize(1,2));

            drawnow
        end

    end
    itt=0;
    o=0;
    while 1==1
        itt=itt+1;
        o=o+1;

        if o>40 && ants > 2
            ants=ants-1;
           
        end
        
         if evaporation <= 0.995
            evaporation = evaporation + 0.005;
            o=0;
         end
    %     ants = round(ants/sqrt(itt)) + 50;
    %     pheromone = itt *pheromone;
    %     
        for ant = 1:ants
        checkList{itt,ant} = [];
        PosList2{itt,ant}=[];

            q=1;
            r=1;

            Pos = fliplr(beginPos);
            PosList{itt,ant}(:,1) = fliplr(beginPos);


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

                    
                    if ((sum(choicesProperty(:,q-1) == 1) + sum(choicesProperty(:,q-1) == 2)) <=2)
                        Maze(Pos(1,1),Pos(1,2))=0;
                    end
                    
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
                
               

                if sum((choicesProperty(:,q-1) == 1) + (choicesProperty(:,q-1) == 2))>3 && q>2 
                    
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

                    probs{itt,ant}(:,q) = (choicesProperty(:,q) == 1);
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
            
            PosList{itt,ant}(:,length(PosList{itt,ant})+1:length(PosList2{itt,ant})+length(PosList{itt,ant})) = PosList2{itt,ant};    
            s=s+1;
            
            if length(PosList{itt,ant}) < length(pathLength)

                pathLength = PosList{itt,ant};
                a = ant;
                b = itt;
                s=0;

            end
            
            if length(PosList{itt,ant})< length(PosBestOwn)
                PosBestOwn = PosList{itt,ant};
                s=0;
            
                
            end


        pathLengthAnt(ant,itt,k)=length(PosList{itt,ant});


        disp(['iteration: ' num2str(itt) ' | ant: ' num2str(ant) ' | # no better solution: ' num2str(s) ' | Best pathlength: ' num2str(length(pathLength)) ' | Current best: ' num2str(length(PosBestOwn)) ' | Current Pathlength: ' num2str(length(PosList{itt,ant}))]);
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
       
        
        if itt>5
            if itt == 100 || s>200 || (length(PosList{itt,1}) == length(PosList{itt-1,1}) && length(PosList{itt,1}) == length(PosList{itt-2,1}) && length(PosList{itt,1}) == length(PosList{itt-3,1})) && length(PosList{itt,1}) == length(PosList{itt-4,1}) && length(PosList{itt,1}) == length(PosList{itt-5,1});
                    
                converged(k) = length(PosList{itt,ant});
                [size1,size2,size3] =size(pathLengthAnt);
                for i=1:size2
                    meanAnt(i)=mean(mean(pathLengthAnt(:,i,:)));
                    drawnow
                end
                    
                break
                
                
                
            end
        end
    end
    
    if k==100;
        break
    end
    

end



% %% Plotting graphical the ammount of pheremones on a certain spot
% 
% % input of the mazes 
% % Maze = dlmread([mazestring ' maze.txt'],'',1,0);
% 
% % input of the starting and ending coordinates of the maze
% Maze_c = load([mazestring ' coordinates.txt']);
% h = figure(2);
% hold on
% set(h, 'Position', [0 0 mazeSize(1,2)*15 mazeSize(1,1)*15]);
% movegui(h,'center')
% 
% for i = 1:mazeSize(1,1)
%     
%     for j = 1:mazeSize(1,2)
%         
%         
%         for k=1:4
%             dir = [i,j] - windDir(k,:);
%             
%             if dir(1,1) ~=0 && dir(1,2) ~=0 && dir(1,1) < mazeSize(1,1) && dir(1,2) < mazeSize(1,2)
%                 
%                 pherclose(k) = pherMatrix(dir(1,1),dir(1,2));
%             
%             else
%                 pherclose(k) = 0;
%             end
%         end
%         
%         pherclose(5) = pherMatrix(i,j);
%         
%         if sum(pherclose == 0) == 5
%         
%             pherclose(1) = inf;
%             
%         end
%         
%         plot(j,-i,'d','Color',[1-pherMatrix(i,j)/max(pherclose),((pherMatrix(i,j)/max(pherclose))),0],'LineWidth',12,'MarkerFaceColor',[1-pherMatrix(i,j)/max(pherclose),((pherMatrix(i,j)/max(pherclose))),0])
%         hold on
%         
%     end
% end
% 
% drawnow
%%
lengte = 0;
lengte = length(pathLength)+lengte-1;

winddir{1}=[num2str(lengte) ';'];
winddir{2}=[num2str(pathLength(2,1)-1) ', ' num2str(pathLength(1,1)-1) ';'];

k=2;
    for i = 1:length(pathLength)-1
          
        winddirr = pathLength(:,i+1)'-pathLength(:,i)';        
        
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

fid = fopen('5_easy.txt','w');

for i=1:2
    fprintf(fid,'%s\n',num2str(winddir{i}));
end

for i=3:length(winddir)
    fprintf(fid,'%s',num2str(winddir{i}));
end
fclose(fid)
