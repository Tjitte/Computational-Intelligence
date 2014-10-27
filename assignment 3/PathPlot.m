function [] = PathPlot(Path,nodes,products)
close all
    %% Plotting the route of the ant

    % input of the mazes 
    Maze = dlmread('Mazes/hard maze.txt','',1,0);
    mazeSize = size(Maze);
    
     h = figure(2);

     hold on

     set(h, 'Position', [0 0 mazeSize(1,2)*15 mazeSize(1,1)*15]);
     movegui(h,'center')

    % setting the axes ranges
    axis([0 81 0 81])
    MatrixPlot = zeros(2,(mazeSize(1,1)+2)*(mazeSize(1,2)+2));

    % looping through the whole maze
    f=0;
    for i=0:mazeSize(1,2)+1;
        for j=0:mazeSize(1,1)+1;

            % if the i or j is not in the matrix create an X
            if i < 1 || i > mazeSize(1,1) || j < 1 || j > mazeSize(1,2)

                f=f+1;
                MatrixPlot(:,f) = [-i+mazeSize(1,2)+1 ; j];

            else

                % if the position on the maze has a wall (0) plot an X there
                if Maze(i,j)==0

                    f=f+1;
                    MatrixPlot(:,f) = [-i+mazeSize(1,2)+1; j];

                end

                if Maze(i,j)==1

                    f=f+1;
                    MatrixPlot(:,f)=[0 ; 0];

                end

            end

        end

    end
        



    h=plot(MatrixPlot(1,1),MatrixPlot(2,1),'xk','linewidth',mazeSize(1,2)/10);
    for f=2:length(MatrixPlot)/(mazeSize(1,1)+2):length(MatrixPlot)
        set(h,'Xdata',MatrixPlot(2,1:f+mazeSize(1,2)),'Ydata',MatrixPlot(1,1:f+mazeSize(1,2)));
        drawnow
    end
    
    for i=1:length(products)
        plot(products(i,1)+1,-products(i,2)+mazeSize(1,2),'dg','MarkerFaceColor','g')
    end

    i=1;
    h=plot(Path{1}(2,i),-Path{1}(1,i)+mazeSize(1,2),'xr','LineWidth',5);
    j=plot(Path{1}(2,i),-Path{1}(1,i)+mazeSize(1,2),'xm','LineWidth',2);
    drawnow
    steps=0;
    for k=1:length(Path)
        for i=2:length(Path{k})
            % plot the route of the ant

            set(h,'Xdata',Path{k}(2,i),'Ydata',-Path{k}(1,i)+mazeSize(1,2)+1);
            set(j,'Xdata',Path{k}(2,1:i),'Ydata',-Path{k}(1,1:i)+mazeSize(1,2)+1);
            
            if length(Path{k}) == i
                plot(Path{k}(2,i),-Path{k}(1,i)+mazeSize(1,2)+1,'dr','MarkerFaceColor','r')
            end
                hold on
                
            steps = steps +1;
            disp(steps);
            drawnow

        end
    end
