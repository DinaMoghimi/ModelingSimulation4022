clc; clear; close all;

%% setup
A=importdata("pHdata.dat");
%data Columns
timeSteps=A(:,1);
u1=A(:,2);  %input
u2=A(:,3);  %input
y=A(:,4);   %output

figure()
scatter(u1,y)
xlabel('u1')
ylabel('y')
grid on
figure()
scatter(u2,y)
xlabel('u2')
ylabel('y')
grid on
figure()
scatter3(u1,u2,y)
xlabel('u1')
ylabel('u2')
zlabel('y')
% normalizing data
u1 = (u1 - min(u1)) / (max(u1) - min(u1));
u2 = (u2 - min(u2)) / (max(u2) - min(u2));
y = (y - min(y)) / (max(y) - min(y));

U=[ones(size(y)),u1,u2];

%% Least square
x = [u1 , u2];
X=[ones(size(y))];
theta=zeros(9,4);
for n=1:4
    X=[X,x.^n];

    for i=1:2*n+1
        LS=lsqr(X,y);
        theta(i,n)=(LS(i,1));
    end
    intercept=theta(1,n);
    slope11=theta(2,n);
    slope21=theta(3,n);
    slope12=theta(4,n);
    slope22=theta(5,n);
    slope13=theta(6,n);
    slope23=theta(7,n);
    slope14=theta(8,n);
    slope24=theta(9,n);
    yFit = intercept + slope11*u1 + slope21*u2 + slope12*u1.^2 + slope22*u2.^2 ...
          + slope13*u1.^3 + slope23*u2.^3 + slope14*u1.^4 + slope24*u2.^4;
    e_LS = y - yFit;
    %plot
    figure(4)
    subplot(2,2,n)
    scatter3(u1,u2,yFit)
    hold on
    scatter3(u1,u2,y,'filled')
    title( num2str(n),'order')
    xlabel('u1')
    ylabel('u2')
    zlabel('y')
    legend('Fitted', 'Data')

    figure(5)
    subplot(2,2,n)
    plot(e_LS)
    title( num2str(n),'order')
    xlabel('Sample Index')
    ylabel('Residual')
    grid on

end

%% forgetting factor
n=0;

for lambda = 0.09:0.05:0.99
    n=n+1;

    theta = zeros(size(U,2),1);
    P = eye(size(U,2)) / lambda;
    
    for i=1:length(y)
        u_i = U(i,:)';
        y_predict = u_i'*theta;
        e = y(i) - y_predict;
        K = P*u_i/(lambda + u_i'*P*u_i);
        theta = theta + K*e;
        P = (P - K*u_i'*P)/lambda;
    end

% error
intercept=theta(1);
slope_u1=theta(2);
slope_u2=theta(3);
yFit_ff = intercept + slope_u1*u1 + slope_u2*u2;
e_ff = y - yFit_ff;
figure(6)
subplot(10,2,n)
plot(e_ff)
grid on
title('\lambda=',num2str(lambda))
end
%% sliding window
window_size=1000;
step_size=1;

num_points=length(y);
num_windows = floor((num_points - window_size) / step_size) + 1;
intercept_sw = zeros(num_windows,1);
slope1_sw = zeros(num_windows,1);
slope2_sw = zeros(num_windows,1);
e_sw = zeros(num_windows,window_size);

for i=1:num_windows
    Start = (i-1) * step_size + 1;
    End = Start + window_size - 1;
    u1_inWindow = u1(Start:End);
    u2_inWindow = u2(Start:End);
    y_inWindow = y(Start:End);

    U_window = [ones(size(y_inWindow)),u1_inWindow,u2_inWindow];
    theta_window = pinv(U_window) * y_inWindow;
    intercept_sw(i)=theta_window(1);
    slope1_sw(i)=theta_window(2);
    slope2_sw(i)=theta_window(3);

    yFit_window = U_window * theta_window;
    e_sw(i,:) = y_inWindow - yFit_window ;
end

%plot
figure(7)
for i=1:num_windows
    Start = (i-1) * step_size + 1;
    End = (i-1) * step_size + window_size ;
    u1_inWindow = u1(Start:End);
    u2_inWindow = u2(Start:End);

    yFit_sw = intercept_sw(i) + slope1_sw(i)*u1_inWindow + slope2_sw(i)*u2_inWindow;
    scatter3(u1,u2,y,'filled','k')
    hold on
    plot3(u1_inWindow,u2_inWindow,yFit_sw)

end
figure(8)
plot(e_sw)
grid on


%% RLS sliding window
window_size=1000;
step_size=1;

num_points=length(y);
intercept_swR = zeros(num_windows,1);
slope1_swR = zeros(num_windows,1);
slope2_swR = zeros(num_windows,1);
e_swR = zeros(num_windows,window_size);

for i = 1 : num_points - window_size + 1
    Start = (i-1) * step_size + 1;
    End = Start + window_size - 1;
    u1_inWindow = u1(Start:End);
    u2_inWindow = u2(Start:End);
    y_inWindow = y(Start:End);

    P = eye(3);
    theta_window = zeros (3,1);
    window_errors = zeros(1,window_size);

    for j =1:window_size
        u=[1; u1_inWindow(j); u2_inWindow(j)];
        e = y_inWindow(j) - u' * theta_window;
        K = (P*u) / (1 + u'*P*u);
        theta_window= theta_window +  K*e ;
        P = P - K*u'*P;
        window_errors(j) = e;
    end

    e_swR(i,:) = window_errors;

end

%plot
figure(9)
for i=1:num_windows
    Start = (i-1) * step_size + 1;
    End = (i-1) * step_size + window_size ;
    u1_Window = u1(Start:End);
    u2_Window = u2(Start:End);
    yFit_swR = intercept_swR(i) + slope1_swR(i)*u1_Window + slope2_swR(i)*u2_Window;
    plot3(u1_Window,u2_Window,yFit_swR)
    hold on
end
scatter3(u1,u2,y,'filled','b')

figure(10)
plot(e_swR)
grid on
