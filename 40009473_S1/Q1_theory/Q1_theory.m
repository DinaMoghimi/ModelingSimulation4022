clc; clear; close all;

% data points
x = [1; 2 ;3]; % input
y = [4; 5 ;7]; % output

X = [ones(size(x)),x]; %  design matrix 

% method 1
coefficients = X \ y;
m = coefficients(2);
b = coefficients(1);
% method 2
theta = (inv(X'*X))*X'*y;
theta_1 = theta(1);
theta_2 = theta(2);
fprintf('theta_{1} = %.2f\n',theta_1)
fprintf('theta_{2} = %.2f\n\n',theta_2)
%% Predict y for a specific x (e.g., x = 4)
x_new = 4 ;
y_predicted = m * x_new + b;

fprintf('Predicted y for x = %.2f is y = %.2f\n', x_new, y_predicted);

% Plot the data points and the regression line
figure;
scatter(x, y, 'o', 'DisplayName', 'Data Points');
hold on;
x_fit = linspace(min(x), max(x), 100);
y_fit = m * x_fit + b;
plot(x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Regression Line');
xlabel('x');
ylabel('y');
title('Least Squares Regression');
legend('Location', 'best');
grid on;

%%  Predict y-values using the model
% Fit a linear regression model
mdl = fitlm(x, y);

% Predict y-values using the model
x_new = 10 * rand(1,100)'; % New x-values for prediction
y_pred = predict(mdl, x_new);

% Plot 
figure()
scatter(x_new,y_pred,'DisplayName', 'Data Points');
hold on;
plot(x_new, y_pred, 'r-', 'LineWidth', 2, 'DisplayName', 'Regression Line');
xlabel('x');
ylabel('y');
title('Least Squares Regression');
legend('Location', 'best');
grid on;

%% 95% Prediction Intervals 
degree = 1; % Degree of the fit
[p,S] = polyfit(x,y,degree);

alpha = 0.05; % Significance level
[yfit,delta] = polyconf(p,x,S,'alpha',alpha);

figure()
plot(x,y,'bo')
hold on
plot(x,yfit,'g-','LineWidth',2)
plot(x,yfit-delta,'r--',x,yfit+delta,'r--')
legend('Data','Fit','95% Prediction Intervals')
xlabel('x');
ylabel('y');
title(['Fit:',texlabel(polystr(round(p,2)))])
hold off