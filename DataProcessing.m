    %% RESET
clear all;
close all;
clc;
format compact;

%% LOAD
Folder = 'ForceData';
Data = {};
Variable = 'Force';

dirList = what(Folder)
for i = 1:1:max(size(dirList.mat))
    File = strcat(Folder, '/', dirList.mat(i));
    Case = load(File{1}, '-mat');
    Data{i} = Case.Polymers;
end

%% PROCESS
L = zeros(size(Data,2), 250);
Weight = zeros(size(Data,2),1);
for i = 1:1:size(Data,2)
    for j = 1:1:size(Data{1,i},2)
        L(i,:) = L(i,:) + sqrt(Data{1,i}{1,j}.BeadPosition(:,1).^2 + Data{1,i}{1,j}.BeadPosition(:,2).^2)';
        Weight(i) = Weight(i) + 10.^Data{1,i}{1,j}.Weight;
    end
    L(i,:) = L(i,:) / size(Data{1,i},2);
    Weight(i) = Weight(i) / size(Data{1,i},2);
end

%% PLOT

h = surf(L);
