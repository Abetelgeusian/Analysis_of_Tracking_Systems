% Harshal Varpe
% ECE 8540
% Lab 4 -  Kalman Filter

clear all
clc
%% Read data files
data_1d = importdata('1D-data.txt');


%%Part 1 - 1D data
T = 1; % time interval of 1 second
phi = [1 T;0 1]; % State transition matrix
Q = [0 0;0 0.0001] ; % Dynamic Noise covariance
R = 1 ; % Measurement noise
M = [1 0]; % Observation Matrix
X_p = [0;0]; % Previous state matrix / initial state matrix
S_p = [1 0;0 1]; % state covariance
plt_data = zeros(1,length(data_1d)); % output data

for t = 1:1:length(data_1d)
    Yt(t) = data_1d(t);
    X_n = phi * X_p ;
    S_n = (phi * S_p * phi') + Q ;
    Kt = S_n * M'/((M*S_n*M')+ R);
    X_p = X_n + Kt * (Yt(t) - M*X_n);
    S_p = (eye(2) - Kt*M) * S_n ;
    plt_data(t) = X_p(1);
end

figure

plot((1:length(data_1d)),data_1d,'k')
hold on;
plot((1:length(data_1d)),plt_data, 'k.-','linewidth',0.5)
axis([100 190 -3 3])
xlabel("Time interval T");
ylabel("Position X");
legend("Measurements","Filtered Estimate");
title("1D Constant Velocity - Ratio Case 2")