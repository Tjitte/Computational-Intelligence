% Function for determining the possible routes for the ant going through a
% maze. ChoicesProperty is a vector of 4 numbers, which contain a 0 1 or a
% 2. the first number is east, the second north, the third west and the
% fourth south.
% 0 - wall
% 1 - undiscovered pathway
% 2 - discovered pathway


function [choicesProperty] = choicesFunc(Pos,Maze,PosList)
                
                choicesProperty = zeros(1,4);
                windDir(1,:)=[0,1];
                windDir(2,:)=[-1,0];
                windDir(3,:)=[0,-1];
                windDir(4,:)=[1,0];
                choices = zeros(4,2);
                
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
                        choicesProperty(i)=Maze(choices(i,1),choices(i,2));

                            % If the choice is 1
                        if choicesProperty(i) == 1;

                           % If so, then change its choiceProperty to 2
                            if sum(PosList{itt,ant}(1, :) == choices(i,1) & PosList{itt,ant}(2, :) == choices(i,2));

                                choicesProperty(i)=2;

                            end

                        end

                    else

                        % if it is not valid, it should be 0
                        choicesProperty(i)=0;

                    end

                end