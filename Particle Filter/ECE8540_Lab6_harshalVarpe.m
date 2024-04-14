clear
clc
close all
% --------------------------------------------------------------------------
data = importdata("magnets_data.txt");

gt_pos = data(:,1);
gt_vel = data(:,2);
gt_sensor = data(:,3);

% Partical Filter Variables
M = 1000; %number of particles
sig_a = 0.0625; % Variance of dynamic noise
sig_m = 4;
sig_n = 0.003906; % variance of measurement noise
Xm1 = -10; % position of magnet
Xm2 = 10;   % position of magnet

% Initialization of weights and state variables
x_pos = zeros(M,1);
x_vel = zeros(M,1);
x_pos_prev = zeros(M,1);
x_vel_prev = zeros(M,1);

wt_norm = 1/M * zeros(M,1);
wt_prev = 1/M * ones(M,1);
wt_up = 1/M * zeros(M,1);

index = zeros(M,1);
x_pos_1 = zeros(M,1);
x_vel_1 = zeros(M,1);
wt_norm_1 = zeros(M,1);

%plot flags
resampling_count = 0
wt_plt = 0; %weight plots
resample_plt = 50; % resampling plots will start at this point

%resampling threshold
rs_thresh = 0.5 ;

X1 = 1:length(gt_sensor);
X2 = 1:M;

result = zeros(length(gt_sensor),1);


q = 2000;
for i = 1:length(gt_sensor)
%     E(i) = 0;

    for j = 1:1:M
        x_pos(j) = x_pos_prev(j) + x_vel_prev(j);
        if(x_pos_prev(j) < -20)
            x_vel(i) = 2;
            
        elseif(x_pos_prev(j) >= -20 && x_pos_prev(j) < 0 )
            x_vel(j) = x_vel_prev(j) + abs(normrnd(0,sig_a));
            
        elseif(x_pos_prev(j) >= 0 && x_pos_prev(j) <= 20)
            x_vel(j) = x_vel_prev(j) - abs(normrnd(0,sig_a));
            
        elseif(x_pos_prev(j) > 20)
            x_vel(j) = -2;
        end
        
        y_t = ( (1/(sqrt(2*pi)*sig_m)) * exp((-(x_pos_prev(j)-Xm1)^2) / (2*sig_m^2)) ) + (1/(sqrt(2*pi)*sig_m)) * exp(-(x_pos_prev(j)-Xm2)^2 / (2 *(sig_m^2)) );
        
        %weight updates via ideal measurement
        p = (1/(sqrt(2*pi)*sig_n)) * exp (-(y_t-gt_sensor(i))^2 / (2 * sig_n^2)) ;
        
        wt_up(j) = wt_prev(j) *  p;        
        wt_prev(j) = wt_up(j);
        x_pos_prev(j) = x_pos(j);
        x_vel_prev(j) = x_vel(j);
    end
    %weight normalization
    wt_norm = wt_up ./ sum(wt_up);
%     sum(wt_up)
    E =0 ;
    for k = 1:M
        E = E + (wt_norm(k) * x_pos(k));
%         E(i) = E(i) + wt_norm(k) + x_pos(k);
    end
    
    result(i) = E;
    
    CV = 0;

    for k = 1:M
        CV = CV + (M*wt_norm(k) - 1)^2;
    end
    
    CV = CV/M;
    ESS = M/(1+CV);
    
    
    if ESS < rs_thresh*M
        
%         if(wt_plt > 0)
%             figure(i)
%             plot(x_pos, wt_norm,'kx')
%             xlabel('Position');
%             ylabel('Weight');
%             title("Before Resampling")
%         end
        
        resampling_count = resampling_count + 1
        if resampling_count >= resample_plt || resampling_count <= resample_plt + 2
            wt_plt = 1;
        end
        
        Q = cumsum(wt_norm);
        t = rand(M+1,1);
        T = sort(t);
        T(M+1) = 1;
        x=1;
        y=1;
        while x<=M
           
            if T(x) < Q(y)
               index(x) = y;
               x = x+1;
           else
               y=y+1;
           end
        end
        for x = 1:M
            x_pos(x)      = x_pos(index(x));
            x_pos_prev(x)  = x_pos_prev(index(x));

            x_vel(x)        = x_vel(index(x));
            x_vel_prev(x)    = x_vel_prev(index(x)) ;

            wt_norm(x)          = 1/M;
            wt_prev(x)     = 1/M;
        end
%            for x=1:M
%                x_pos_1(x) = x_pos(index(x));
%                wt_norm_1(x) = 1/M;
%                x_vel_1(x) = x_vel(index(x));
%            end
%         x_pos = x_pos_1;
%         wt_norm = wt_norm_1;
%         x_vel = x_vel_1;
%         if(wt_plt > 0)
%             figure(i+100)
%             plot(x_pos, wt_norm,'kx')
%             xlabel('Position');
%             ylabel('Weight');
%             title("After Resampling")
% %             q = q+1;
%         end
    end
%     wt_prev = wt_norm; % uncmnt this blokc it works
%     x_pos_prev = x_pos;
%     x_vel_prev = x_vel;

    if resampling_count > resample_plt
        wt_plt = 0;
    end


end
% figure(i+1)
% plot(gt_pos,"kx")
% hold on
% plot(result,"k-")
% xlabel("Time")
% ylabel("Position")
% legend("Ground truth","Estimated Position")

hold off
figure(i+2)
plot(gt_sensor,"k")
xlabel("Time")
ylabel("Field Strength")