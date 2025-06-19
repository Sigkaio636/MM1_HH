%% Test the reduction hypothesis for model 2
close all
%% Check bistability
figure(1); hold on; grid on;
tspan = [0 50];
for midx = 0.1:0.2:0.9
    for hidx = 0.1:0.2:0.9
        for nidx = 0.1:0.2:0.9
            for vidx = 0:1
                x0idx = [vidx; nidx; midx; hidx; 10];
                [tidx, xidx] = ode45(@HHfull1, tspan, x0idx);
                plot(tidx, xidx(:,1))
            end
        end
    end
end
%% Compare trajectories 
x0 = [0.05; 0.32; 0.05; 0.59; 4]; % Initial state <- Play changing the intensity
tspan = [0 25];
[tf, xf] = ode45(@HHfull2, tspan, x0); % Solve for trajectory for full model
%xr0 = [x0(1); x0(2); x0(5)]; 
%[tr, xr] = ode45(@HHredu1, tspan, xr0); % Solve for trajectory for reduced model

% Compare trajectories graphically
figure(2)
plot(tf, xf(:,1), 'b-', tr, xr(:,1), 'b--');
xlabel('Time');
ylabel('Action Potential');

%am = @(V) 0.08 .* (V + 56) ./ (1 - exp(-(V + 56) / 6.8));
%bm = @(V) 0.8 .* exp(-(V + 56) / 18);
%minf =@(V) am(V) ./ (am(V) + bm(V));

%figure(2)
%hold on; 
%plot(tf, xf(:,2), 'r-', tf, xf(:,3), 'g-', tf, xf(:,4), 'b-');
%plot(tr, xr(:,2), 'r--', tr, minf(xr(:,1)), 'g--', tr, (0.89 - 1.1 * xr(:,2)), 'b--');
%xlabel('Time');
%ylabel('Gate activation');

% We can see that qualitatively they behave quite similar
