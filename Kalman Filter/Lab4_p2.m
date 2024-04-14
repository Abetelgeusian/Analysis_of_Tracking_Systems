% Harshal Varpe
% ECE 8540
% Lab 4 -  Kalman Filter

clear all
clc
%% Read data files
data_2d = importdata('2D-UWB-data.txt');

%%Part 1 - 1D data
T = 1; % time interval of 1 second
phi = [1 0 T 0;0 1 0 T; 0 0 1 0;0 0 0 1]; % State transition matrix

Q = [0 0 0 0;0 0 0 0;
    0 0 0.1 1; 0 0 1 0.1] ; % Dynamic Noise covariance

R = [10 0; 0 10] ; % Measurement noise
M = [1 0 0 0; 0 1 0 0]; % Observation Matrix

X_p = [data_2d(1,1);data_2d(1,2);0;0]; % Previous state matrix / initial state matrix
S_p = eye(4); % state covariance
plt_data = zeros(length(data_2d(:,1)),2); % output data

for t = 1:1:length(data_2d(:,1))
    Yt(1,1) = data_2d(t,1);
    Yt(2,1) = data_2d(t,2);
    X_n = phi * X_p ;
    S_n = (phi * S_p * phi') + Q ;
    Kt = S_n * M'/((M*S_n*M')+ R);
    X_p = X_n + Kt * (Yt - M*X_n);
    S_p = (eye(4) - Kt*M) * S_n ;
    plt_data(t,1) = X_p(1,1);
    plt_data(t,2) = X_p(2,1);
end
% 
figure
plot((1:length(data_2d(:,1))),data_2d(:,1),'k')
hold on;
plot((1:length(data_2d)),plt_data(:,1) , 'k.-','linewidth',0.5)
% axis([20 50 270 320])
xlabel("Time interval T");
ylabel("Position X");
legend("Measurements","Filtered Estimate");
title("2D Constant Velocity")


figure
plot((1:length(data_2d(:,1))),data_2d(:,2),'k')
hold on;
plot((1:length(data_2d)),plt_data(:,2) , 'k.-','linewidth',0.5)
% axis([10 30 600 650])
xlabel("Time interval T");
ylabel("Position Y");
legend("Measurements","Filtered Estimate");
title("2D Constant Velocity")




