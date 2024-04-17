


clc;clear;close all


file_path='batch-yield-and-purity.dat';
data_table=readtable(file_path);

x = data_table(:,1).Variables;
y = data_table(:,2).Variables;

% Normalize x1 and y
x = (x - min(x)) / (max(x) - min(x));
y = (y - min(y)) / (max(y) - min(y));

% Scatter plot of the original data
figure(1)
scatter(x,y);
xlabel('x');
ylabel('y');
title('Scatter Plot: Original Data');
grid on;
% Linear regression
theta = pinv(x) * y;

% Calculate the fitted values
y_fit = theta * x;
% Scatter plot of the original data and the fitted values
figure;
scatter(x, y, 'filled')
hold on
scatter(x, y_fit, 'r', 'filled');
xlabel('x1');
ylabel('y');
title('Scatter Plot: Original Data vs Fitted Values');
grid on;
legend('Original Data', 'Fitted Values');
% Calculate the error
E = y - y_fit;

% Error plot
figure;
plot(E, 'k', 'LineWidth', 1.5);
xlabel('Data Point Index');
ylabel('Error (y - y_{fit})');
title('Error Plot: Actual Output - Fitted Output');
grid on;