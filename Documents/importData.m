function [TrainingSamplesOfFeatures, TrainingTargetsOfSamples] = importData()

% import features for training
TrainingSamplesOfFeatures = csvread(strcat(pwd,'/Documents/Data/features.txt'));

% import targets for training
TrainingTargetsOfSamples = csvread(strcat(pwd,'/Documents/Data/targets.txt'));

% import features to calculate targets
% SamplesOfFeatures = csvread(strcat(pwd,'\Data\unknown.txt'));