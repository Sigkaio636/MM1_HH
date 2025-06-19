%% Test the reduction hypothesis for model 1
close all
%% Compare trajectories 
x0 = [0.05; 0.32; 0.05; 0.59; 10]; % Initial state <- Play changing the intensity
tspan = [0 50];
[tf, xf] = ode45(@HHfull1, tspan, x0); % Solve for trajectory for full model
xr0 = [x0(1); x0(2); x0(5)]; 
[tr, xr] = ode45(@HHredu1, tspan, xr0); % Solve for trajectory for reduced model

% Compare trajectories graphically
figure(1)
plot(tf, xf(:,1), 'b-', tr, xr(:,1), 'b--');
xlabel('Time t ');
ylabel('Action Potential V');

am = @(V) 0.1 * (25 - V) ./ (exp(2.5 - V / 10) - 1);
bm = @(V) 4 * exp(-V / 18);
minf =@(V) am(V) ./ (am(V) + bm(V));

figure(2)
hold on; 
plot(tf, xf(:,2), 'r-', tf, xf(:,3), 'g-', tf, xf(:,4), 'b-');
plot(tr, xr(:,2), 'r--', tr, minf(xr(:,1)), 'g--', tr, (0.8882 - 1.041 * xr(:,2)), 'b--');
xlabel('Time t');
ylabel('Gate activation');

% We can see that qualitatively they behave quite similar

%% Find n-h linear regression
Iarr = -5:1:20;
linreg_data = [];
for Iidx = Iarr % iterate for Intensity value
    n_h_data = [];
    for midx = 0.1:0.2:0.9
        for vidx = -10:20:50
            for nidx = 0.1:0.1:0.9
                for hidx = 0.1:0.1:0.9 % iterate for several initial states
                    xidx = [vidx; nidx; midx; hidx; Iidx];
                    [~, x] = ode45(@HHfull1, tspan, xidx);
                    n_h_data = [n_h_data; [x(300:end,2), x(300:end,4)]];
                        % The 300 is here to exclude the transitory points 
                        %  as some initial states are unrepresentative.
                        %  In figure(3) we can visualize this.
                    % solve the trajectory and append the n- h- arrays
                end 
            end
        end 
    end 
    
    % As we considered soo many initial states, apply cross validation by
    %  repetition to estimate the distribution of regression parameters.
    Nrow = size(n_h_data, 1);
    reg_slope_data = [];
    reg_inter_data = [];
    for i = 1:500
        Nsel = round(0.1 * Nrow); % select randomply 10% of data.
        selidx = randperm(Nrow, Nsel);
        n_h_scatter = n_h_data(selidx, :);
        n_h_scatter = n_h_scatter(~any(isnan(n_h_scatter), 2), :);
        % for some reason, Nan values are created. But they are minority 
        %  thus we don't select them
        
        n_tofit = n_h_scatter(:,1);
        h_tofit = n_h_scatter(:,2);
        mean_n = mean(n_tofit);
        mean_h = mean(h_tofit);
        covM_nh = cov(n_tofit, h_tofit);
        reg_slope_data = [reg_slope_data, covM_nh(1,2)/covM_nh(1,1)];
        reg_inter_data = [reg_inter_data, mean_h-mean_n*covM_nh(1,2)/covM_nh(1,1)];
        % calculate the slope and intercept of this subset and append it 
    end 

   linreg_data = [linreg_data, [mean(reg_slope_data), var(reg_slope_data) ; mean(reg_inter_data), var(reg_inter_data)] ];
    % store the data of regression parameters as mean and variance 
    %  for the given Intensity value
end

figure(3)
% Last n_tofit & h_tofit, make a scatterplot to visualize 
hold on;
scatter(n_tofit, h_tofit)
plot(0.1:0.05:0.9, ((0.1:0.05:0.9)-mean_n)*reg_slope_data(end) + mean_h)

% Histograms of regression data for last Intensity value
figure(4)
histogram(reg_slope_data, 30);
figure(5)
histogram(reg_inter_data, 30);
% Clearly the distributions are quite normal

% Calculate the mean value across all Intensity values
estim1_slope = mean(linreg_data(1,1:2:end));
estim1_slope_std = sqrt(sum(linreg_data(1,2:2:end)))./length(linreg_data(1,2:2:end));
estim1_inter = mean(linreg_data(2,1:2:end));
estim1_inter_std = sqrt(sum(linreg_data(2,2:2:end)))./length(linreg_data(2,2:2:end));

% Plot the mean with error bars for the slope and the intercept
figure(6)
bar(Iarr, linreg_data(1,1:2:end)); hold on;
errorbar(Iarr, linreg_data(1,1:2:end), sqrt(linreg_data(1,2:2:end)), '.k', 'LineWidth', 1.5);
yline(estim1_slope, 'r-');
yline(estim1_slope+estim1_slope_std, 'r--');
yline(estim1_slope-estim1_slope_std, 'r--');

figure(7)
bar(Iarr, linreg_data(2,1:2:end)); hold on;
errorbar(Iarr, linreg_data(2,1:2:end), sqrt(linreg_data(2,2:2:end)), '.k', 'LineWidth', 1.5);
yline(estim1_inter, 'r-');
yline(estim1_inter+estim1_inter_std, 'r--');
yline(estim1_inter-estim1_inter_std, 'r--');

