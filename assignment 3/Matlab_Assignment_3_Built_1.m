%% Variables

clear all
close all 
clc

%% Input Maze

% input of the mazes 
easyMaze = dlmread('Mazes\easy maze.txt','',1,0);
mediumMaze = dlmread('Mazes\medium maze.txt','',1,0);
hardMaze = dlmread('Mazes\hard maze.txt','',1,0);

% input of the starting and ending coordinates of the maze
easyMaze_c = load('Mazes\easy coordinates.txt',',');
mediumMaze_c = load('Mazes\medium coordinates.txt',',');
hardMaze_c = load('Mazes\hard coordinates.txt',',');

