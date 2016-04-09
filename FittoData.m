    %% RESET
clear all;
close all;
clc;
format compact;

%% LOAD
Data = load('Data.mat');

R = Data.R;
F = Data.F;
Data = [F',R'];
Data = sortrows(Data);
F = Data(3:end,1);
R = Data(3:end,2);
f = @(c) log10(250 - (250-250.^0.75).* exp(-F.^0.1.*c(1) - F.^0.5))
fc = @(c) sum((R - f(c)).^2);
c0 = [0.3];
options = optimset('MaxFunEvals', 1E8, 'MaxIter', 1E4, 'TolX',1e-12, 'TolFun', 1e-12);
[Y,V] = fminsearch(fc, c0, options)

figure();
hold on;
plot(F,R, '.k');
plot(F, f(Y), 'xr');