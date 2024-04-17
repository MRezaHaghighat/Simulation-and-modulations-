

clc;clear;close all


file_path='pHdata.dat';
data_table=readtable(file_path);

t=data_table(:,1).Variables;
x1 = data_table(:,2).Variables;
x2 = data_table(:,3).Variables;
y = data_table(:,4).Variables;

x1 = (x1 - min(x1)) / (max(x1) - min(x1));
x2 = (x2 - min(x2)) / (max(x2) - min(x2));
y = (y - min(y)) / (max(y) - min(y));

X =[ones(size(x1)).*x1 x2];
% Scatter plot of the original data
figure(1)
scatter(x1,y);
figure(2)
scatter(x2,y)
figure(3)
scatter3(x1,x2,y,'filled');
grid on;
% Linear regression
theta = pinv(X) * y;
thetal=lsqr(X,y);

slop_x1=theta(1);
slop_x2=theta(2);

% Calculate the fitted values
y_fit =slop_x1*x1+slop_x2*x2;
% Scatter plot of the original data and the fitted values
figure;
scatter3(x1,x2, y, 'filled')
hold on
scatter3(x1,x2, y_fit, 'r', 'filled');
xlabel('Acid and Base solution flow in liters');
ylabel('pH of the solution in the tank');
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

%% forgetting factor


lambda = 0.9;  % recommended: 0.7 < lambda < 0.9

% Forgetting Factor
theta = zeros(size(X, 2), 1);
P = eye(size(X, 2)) / lambda;

% Loop through all available data points (up to the length of y)
for i = 1:length(y)
    x_i = X(i, :)';
    y_predicted = x_i' * theta;
    e = y(i) - y_predicted;
    K = P * x_i / (lambda + x_i' * P * x_i);
    theta = theta + K * e;
    P = (P - K * x_i' * P) / lambda;
end

% Error
% Calculate prediction values
y_fit = X * theta;
E = y - y_fit;

% Plotting Error
figure;
plot(E, 'k', 'LineWidth', 1.5);
xlabel('Data point Error');
ylabel('Index');
title('Error Plot = Actual Value - Fitted Value with Forgetting Factor');
grid on;



%% sliding window


% Sliding Window
Window_Size = 50;
Step_Size = 1;

num_points = length(y);
num_windows = (num_points - Window_Size) / Step_Size + 1;
intercepts = zeros(num_windows, 1);
slope_x1 = zeros(num_windows, 1);
slope_x2 = zeros(num_windows, 1);

errors = zeros(num_windows, Window_Size);

% Sliding window least square
for i = 1:num_windows
    start_idx = (i - 1) * Step_Size + 1;
    end_idx = start_idx + Window_Size - 1;

    x1_window = x1(start_idx:end_idx);
    x2_window = x2(start_idx:end_idx);
    y_window = y(start_idx:end_idx);

    X_window = [ones(size(x1_window)), x1_window, x2_window];
    theta_window = pinv(X_window) * y_window;
    intercepts(i) = theta_window(1);
    slope_x1(i) = theta_window(2);
    slope_x2(i) = theta_window(3);

    y_fit_window = X_window * theta_window;

    errors(i, :) = y_window - y_fit_window;
end

% Plotting
figure;
ax1 = subplot(2, 1, 1);
scatter3(x1, x2, y, 'blue', 'filled');
hold on;
for i = 1:num_windows
    x1_window = x1((i - 1) * Step_Size + 1:(i - 1) * Step_Size + Window_Size);
    x2_window = x2((i - 1) * Step_Size + 1:(i - 1) * Step_Size + Window_Size);
    y_fit = intercepts(i) + slope_x1(i) * x1_window + slope_x2(i) * x2_window;
    plot3(x1_window, x2_window, y_fit, 'red');
end
hold off;
title('Fitted Lines using Sliding Window Least Squares Error');
xlabel('Acid solution flow in liters');
ylabel('Base solution flow in liters');
zlabel('pH of the solution in the tank');

ax2 = subplot(2, 1, 2);
plot(errors', 'k', 'LineWidth', 1.5);
xlabel('Data point index');
ylabel('Error (y - y-fit)');
title('Error (y - y-fit)');
grid on;

linkaxes([ax1, ax2], 'x');


%% RLS Method in sliding window

Window_Size = 50;
Step_Size = 1;

% Parameters
intercepts = zeros(length(y) - Window_Size + 1, 1);
slope_x1 = zeros(length(y) - Window_Size + 1, 1);
slope_x2 = zeros(length(y) - Window_Size + 1, 1);

errors = zeros(length(y) - Window_Size + 1, Window_Size);

for i = 1:length(y) - Window_Size + 1
    x1_window = x1(i:i + Window_Size - 1);
    x2_window = x2(i:i + Window_Size - 1);
    y_window = y(i:i + Window_Size - 1);

    P = eye(3);
    theta = zeros(3, 1);
    window_errors = zeros(Window_Size, 1);

    for j = 1:Window_Size
        x = [1; x1_window(j); x2_window(j)];
        e = y_window(j) - x' * theta;
        K = P * x / (1 + x' * P * x);
        theta = theta + K * e;
        P = P - K * x' * P;
        window_errors(j) = e;
    end

    errors(i, :) = window_errors';
    intercepts(i) = theta(1);
    slope_x1(i) = theta(2);
    slope_x2(i) = theta(3);
end

% Plotting
figure;
ax1 = subplot(2, 1, 1);
scatter3(x1, x2, y, 'blue', 'filled');
hold on;
for i = 1:length(y) - Window_Size + 1
    x1_window = x1(i:i + Window_Size - 1);
    x2_window = x2(i:i + Window_Size - 1);
    y_fit = intercepts(i) + slope_x1(i) * x1_window + slope_x2(i) * x2_window;
    plot3(x1_window, x2_window, y_fit, 'red');
end
hold off;
title('Fitted Lines using Sliding Window RLS Method');
xlabel('Acid solution flow in liters');
ylabel('Base solution flow in liters');
zlabel('pH of the solution in the tank');

ax2 = subplot(2, 1, 2);
plot(errors', 'k', 'LineWidth', 1.5);
xlabel('Data point index');
ylabel('Error (y - y-fit)');
title('Errors plot: Actual Value - Fitted Value (Sliding Window RLS Method)');
grid on;

linkaxes([ax1, ax2], 'x');

