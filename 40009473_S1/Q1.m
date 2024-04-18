clc; clear; close all;

%% import data
A=readmatrix('batch-yield-and-purity.csv');

x=A(:,2);
y=A(:,1);
X=[ones(size(x)),x];

figure()
scatter(x,y)
title('data points')
xlabel('x');
ylabel('y');
grid on
hold on
%% 
% without commands
Theta_ncmd1=(inv(X'*X))*X'*y
Theta_ncmd2=pinv(X)*y
% with command
Theta_cmd=lsqr(X,y)

% curve fitting
ytilda = X*Theta_cmd ;
plot(x,ytilda,'r',LineWidth=3)
xlabel('x');
ylabel('y');
title('curve fitting')

% error
error = y - ytilda;
J = error' * error;

%% BLUE

theta = X \ y;
theta_0 = theta(1);
theta_1 = theta(2)

% Display the estimated linear model
fprintf('Estimated linear model: y = %.4f + %.4fx\n', theta_0, theta_1);

% Plot the data points and the fitted line
figure;
scatter(x, y, 'o', 'DisplayName', 'Data Points');
hold on;
x_fit = linspace(min(x), max(x), 100);
y_fit = theta_0 + theta_1 * x_fit;
plot(x_fit, y_fit, 'r-', 'LineWidth',2, 'DisplayName', 'Fitted Line');
xlabel('x');
ylabel('y');
title('Best Linear Unbiased Estimation (BLUE)');
legend('Location', 'best');
grid on;

