%% RESET
clear all;
close all;
clc;
format compact;

%% LOAD
Data = load('Extra/Survivors_ForcePolymer_3.mat');
Polymers = Data.Polymers;
figure();
hold on;
for i = 1:1:size(Polymers,2)
    plot(Polymers{i}.BeadPosition(:,1), Polymers{i}.BeadPosition(:,2), '-o');
end