function [Samples] = importDataGC()

% import features to calculate targets
Samples = csvread(strcat(pwd,'\Data\unknown.txt'));