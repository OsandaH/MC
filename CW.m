% Monte Carlo integration for I = ∫∫ (x^2 + y) dy dx
clear; clc; rng('shuffle');

% Integration bounds
x_min = 2; 
x_max = 3;

N = 1e6;   % number of random points
x_rand = x_min + (x_max - x_min) * rand(N,1);   % uniform in [2,3]
y_rand = 2*x_min^3 + (2*x_max^3 - x_min) * rand(N,1); % trial range for y

% Define region: y between [x , 2x^3]
inside = (y_rand >= x_rand) & (y_rand <= 2*x_rand.^3);

% Evaluate integrand (only for points inside region)
f_vals = (x_rand.^2 + y_rand) .* inside;

% Bounding box area
Area_box = (x_max - x_min) * (2*x_max^3 - x_min);

% Monte Carlo estimate
I_MC = Area_box * mean(f_vals);

fprintf('Monte Carlo estimate: %.2f\n', I_MC);
fprintf('Exact value: 790.55\n');
fprintf('Error: %.2f %%\n', 100*abs(I_MC-790.55)/790.55);

% --------- Plot the region ---------
figure; hold on;
fplot(@(x) x, [x_min x_max], 'r', 'LineWidth',2);         % lower bound y = x
fplot(@(x) 2*x.^3, [x_min x_max], 'b', 'LineWidth',2);    % upper bound y = 2x^3
xlabel('x'); ylabel('y');
title('Integration Region');
legend('y = x','y = 2x^3','Location','northwest');



%%%%%%%% Particle exchange simulation
clear; clc; 
rng('shuffle');

N = 100;          % total number of particles
nl = N;           % initial number on left side
nr = 0;           % initial number on right side
T = 500;          % number of time steps

nl_record = zeros(1,T); % to record particles on left
nr_record = zeros(1,T); % to record particles on right

for t = 1:T
    r = rand;  % random number in [0,1]
    if r < nl/N
        % particle moves left -> right
        nl = nl - 1;
        nr = nr + 1;
    else
        % particle moves right -> left
        nl = nl + 1;
        nr = nr - 1;
    end
    
    % record values
    nl_record(t) = nl;
    nr_record(t) = nr;
end

% Plot results
figure;
plot(1:T, nl_record, 'b', 'LineWidth',1.5); hold on;
plot(1:T, nr_record, 'r', 'LineWidth',1.5);
xlabel('Time step'); ylabel('Number of particles');
legend('Left side','Right side');
title('Particle exchange simulation');




%%%%%%%%%%%%%% Sampling from f(theta) = 1 / (sin^2(theta) + a*cos^2(theta))

clear; clc; rng('shuffle');

N = 10000;               % number of samples
theta_range = [0, 2*pi]; % domain

% Try two values of 'a'
a_values = [0.5, 0.001];

for k = 1:length(a_values)
    a = a_values(k);

    % Define PDF (unnormalized)
    f = @(theta) 1 ./ (sin(theta).^2 + a*cos(theta).^2);

    % Find maximum value of f(theta) for rejection sampling
    theta_grid = linspace(0, 2*pi, 10000);
    f_max = max(f(theta_grid));

    % Generate samples using rejection method
    samples = zeros(N,1);
    count = 0;
    while count < N
        theta_try = 2*pi*rand;            % uniform candidate
        y = f_max * rand;                 % uniform in [0, f_max]
        if y < f(theta_try)
            count = count + 1;
            samples(count) = theta_try;
        end
    end

    % -------- Plot results --------
    figure;
    histogram(samples, 50, 'Normalization','pdf'); hold on;

    % Normalize f(theta) to make it a proper PDF
    Z = trapz(theta_grid, f(theta_grid)); % numerical normalization constant
    f_norm = f(theta_grid)/Z;

    plot(theta_grid, f_norm, 'r', 'LineWidth',2);
    xlabel('\theta'); ylabel('PDF');
    title(sprintf('Sampling from f(\\theta), a = %.3f', a));
    legend('Histogram of samples','Normalized f(\theta)');
end
