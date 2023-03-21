% Housekeeping
clear
clc
close all

% Load good parameters
load goodParams.mat

% Modify data to be useful
trunc = goodParams(:,[5,7,8,9]);

% Plot
figure()
scatter(trunc(:,1)+0.08,trunc(:,2), "filled", "sk", 'DisplayName', 'Top cover mass')
xlabel('Root rib spacing (m)',"Interpreter","latex", "FontSize", 16)
ylabel('Mass (Kg)',"Interpreter","latex", "FontSize", 16)
hold on
scatter(trunc(:,1)+0.08,trunc(:,3), "filled", "^r", 'DisplayName','Rib mass')
scatter(trunc(:,1)+0.08,trunc(:,4), "filled", "og", 'DisplayName', 'Combined mass')
hold off
legend("Interpreter","latex", "FontSize", 16)
grid on
grid minor
set(gca, "FontSize", 16, "TickLabelInterpreter", "latex")