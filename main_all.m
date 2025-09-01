%% Q1.1
clear all; clc;

% Bounding box
xmin = -10; xmax = 10;
ymin = -10; ymax = 10;
A_box = (xmax-xmin)*(ymax-ymin);

N_values = [1000, 10000, 100000];
area_est = zeros(size(N_values));
error_est = zeros(size(N_values));

for k=1:length(N_values)
    N = N_values(k);
    x = xmin + (xmax-xmin)*rand(1,N);
    y = ymin + (ymax-ymin)*rand(1,N);

    % Inequality test
    inside = ( (x.^2 + y.^2 - 2.*x).^2 <= 4*(x.^2 + 4*y.^2) );

    % Fraction of points inside
    p = sum(inside)/N;

    % Estimated area
    area_est(k) = p*A_box;

    % Error estimate (binomial error scaled by box area)
    error_est(k) = A_box*sqrt(p*(1-p)/N);
end

disp('N, Area Estimate, Error');
disp([N_values' area_est' error_est']);

inside = ( (x.^2 + y.^2 - 2.*x).^2 <= 4*(x.^2 + 4*y.^2) );

figure;
scatter(x(inside), y(inside), 1, 'b'); hold on;
scatter(x(~inside), y(~inside), 1, 'r');
xlabel('x'); ylabel('y');
title('Monte Carlo Sampling of Region');
axis equal;

%% Q1.2

% Region: %
% -5/4 - 3/(3+x^2) <= y <= 5/4 + 3/(3+x^2)
% -1/3 - 5*y^2/6 + y^4/6 <= x <= 11/3 - (2/3)sqrt(2)|y|^(3/2)

clc
clear all;

N_val = [100,1000,10000];

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
box_area = (xmax-xmin)*(ymax-ymin);

for N = N_val
    % Generate random points
    x = xmin + (xmax-xmin)*rand(N,1);
    y = ymin + (ymax-ymin)*rand(N,1);

    y_lower = -5/4 -(3./(3+x.^2));
    y_upper = 5/4 + (3./(3+x.^2));

    x_lower = -1/3 - (5*y.^2)/6 + (y.^4)/6;
    x_upper = 11/3 - (2/3)*sqrt(2)*abs(y).^(3/2);
    
    inside = (y_lower <= y & y <= y_upper & x_lower <=x & x<=x_upper);

    area = box_area*sum(inside)/N;

    p = sum(inside)/N;
    error = area*sqrt(p*(1-p)/N);

    fprintf('N=%d, Estimated Area = %.4f, Error ≈ %.4f\n', N, area, error);

    figure;
    scatter(x(inside), y(inside), 5, 'b'); hold on;
    scatter(x(~inside), y(~inside), 5, 'r');
    title('Region (i): Monte Carlo Points');
    xlabel('x'); ylabel('y'); axis equal;

end

%% Q2 MC

clc; 
clear; 
close all;

% Function to integrate
f = @(x) exp(-x.^2);

% True value for reference
true_val = integral(f, 0, 1);

% Number of samples
N = 1000; % change to test different sizes

% Random uniform samples
x_rand = rand(N,1);
I_mc = mean(f(x_rand));

% Error
err_mc = abs(I_mc - true_val);

fprintf('Plain Monte Carlo:\n');
fprintf('Estimate = %.8f, Error = %.2e\n', I_mc, err_mc);

%% Q2 SC 

clc; clear; close all;

f = @(x) exp(-x.^2);
true_val = integral(f, 0, 1);

N = 1000; % number of strata

edges = linspace(0, 1, N+1);
strat_samples = zeros(N,1);

for i = 1:N
    u = rand;
    x_strat = edges(i) + u*(edges(i+1) - edges(i));
    strat_samples(i) = f(x_strat);
end

I_strat = mean(strat_samples);
err_strat = abs(I_strat - true_val);

fprintf('Stratified Sampling:\n');
fprintf('Estimate = %.8f, Error = %.2e\n', I_strat, err_strat);

%% Q2.adaptive

clc; clear;

N = 5000;              % total number of samples
M = 10;                % number of strata (bins)
true_val = integral(@(x) exp(-x.^2), 0, 1); % reference value

% Step 1: define bin edges
edges = linspace(0,1,M+1);

% Step 2: pilot sampling (say 20 samples/bin)
pilotN = 20;
var_est = zeros(M,1);

for j = 1:M
    a = edges(j); b = edges(j+1);
    x_pilot = a + (b-a)*rand(pilotN,1);
    f_pilot = exp(-x_pilot.^2);
    var_est(j) = var(f_pilot);  % variance estimate
end

% Step 3: allocate samples proportional to variance
alloc = round((var_est/sum(var_est)) * (N - M*pilotN));
alloc = max(alloc,1); % at least 1 per bin

% Step 4: final sampling
I_adaptive = 0;
for j = 1:M
    a = edges(j); b = edges(j+1);
    Nj = alloc(j);
    x_samples = a + (b-a)*rand(Nj,1);
    f_samples = exp(-x_samples.^2);
    I_adaptive = I_adaptive + mean(f_samples) * (b-a);
end

fprintf('Adaptive Sampling Result = %.6f | True Value = %.6f | Error = %.6f\n',...
    I_adaptive, true_val, abs(I_adaptive-true_val));

%% Q2.Impotance

clc; clear;

N = 5000;               % number of samples
true_val = integral(@(x) exp(-x.^2), 0, 1); % reference value

% ---------- Importance Sampling ----------
% Choose weight function w(x) = 2*(1-x), 0<=x<=1
w = @(x) 2*(1-x);

% CDF of w(x): CDF(x) = 2*(x - x^2/2)  => invert numerically
% Solve y = 2*(x - x^2/2) for x using quadratic formula
y = rand(N,1);  % uniform random numbers in [0,1]
x_samples = 1 - sqrt(1 - y); % inverse CDF method

% Evaluate f(x)/w(x)
f_samples = exp(-x_samples.^2);
I_importance = mean(f_samples ./ w(x_samples));

fprintf('Importance Sampling Result = %.6f | True Value = %.6f | Error = %.6f\n',...
    I_importance, true_val, abs(I_importance-true_val));

%% Q4

% Breit-Wigner Distribution Event Generator - Script Version
% Compatible with MATLAB R2013b

clear all;
close all;
clc;

% Parameters
M = 748; % Mass in MeV (m in the formula)
Gamma = 12; % Width in MeV (? in the formula)
N_events = 1000; % Number of events to generate

% Mass range for generation (extend beyond the peak)
m_min = M - 5*Gamma;
m_max = M + 5*Gamma;

fprintf('Breit-Wigner Distribution Event Generator\n');
fprintf('=========================================\n');
fprintf('Peak mass (M): %.1f MeV\n', M);
fprintf('Width (Gamma): %.1f MeV\n', Gamma);
fprintf('Events to generate: %d\n\n', N_events);

%% Part (i): Generate 1000 events uniformly using weight factor B
fprintf('Part (i): Generating events using weight factor method\n');
fprintf('-----------------------------------------------------\n');

% Initialize arrays
accepted_masses_weighted = [];
attempt_count = 0;

% Calculate B_max for normalization
m_test_range = linspace(m_min, m_max, 10000);
Gamma_half = Gamma/2;
B_test_range = Gamma_half ./ ((m_test_range - M).^2 + Gamma_half^2);
B_max = max(B_test_range);

% Accept/reject method using weights
while length(accepted_masses_weighted) < N_events
    % Generate random mass
    m_test = m_min + (m_max - m_min) * rand();
    
    % Calculate Breit-Wigner weight factor B
    B_test = Gamma_half / ((m_test - M)^2 + Gamma_half^2);
    
    % Accept with probability B_test/B_max
    if rand() < B_test/B_max
        accepted_masses_weighted = [accepted_masses_weighted; m_test];
    end
    attempt_count = attempt_count + 1;
end

fprintf('Generated %d events using weight factor method\n', length(accepted_masses_weighted));
fprintf('Efficiency: %.2f%% (%d attempts for %d events)\n\n', ...
    100*N_events/attempt_count, attempt_count, N_events);

%% Part (ii): Accept/reject method with B/B? condition
fprintf('Part (ii): Accept/reject method with B/B? > random condition\n');
fprintf('----------------------------------------------------------\n');

% Calculate B? (value of B when m=M)
B_0 = Gamma_half / ((M - M)^2 + Gamma_half^2);
B_0 = Gamma_half / Gamma_half^2;
B_0 = 1/Gamma_half;

% Initialize arrays for second method
accepted_masses_reject = [];
attempt_count_2 = 0;

while length(accepted_masses_reject) < N_events
    % Generate random mass
    m_test = m_min + (m_max - m_min) * rand();
    
    % Calculate Breit-Wigner weight factor B
    B_test = Gamma_half / ((m_test - M)^2 + Gamma_half^2);
    
    % Accept if B/B? > random number
    if (B_test/B_0) > rand()
        accepted_masses_reject = [accepted_masses_reject; m_test];
    end
    attempt_count_2 = attempt_count_2 + 1;
end

fprintf('Generated %d events using accept/reject method\n', length(accepted_masses_reject));
fprintf('Efficiency: %.2f%% (%d attempts for %d events)\n\n', ...
    100*N_events/attempt_count_2, attempt_count_2, N_events);

%% Plot results
figure('Position', [100, 100, 1200, 400]);

% Plot part (i) results
subplot(1, 3, 1);
hist(accepted_masses_weighted, 50);
title('Part (i): Weight Factor Method');
xlabel('Mass (MeV)');
ylabel('Frequency');
grid on;

% Plot part (ii) results
subplot(1, 3, 2);
hist(accepted_masses_reject, 50);
title('Part (ii): Accept/Reject Method');
xlabel('Mass (MeV)');
ylabel('Frequency');
grid on;

% Plot theoretical Breit-Wigner curve
subplot(1, 3, 3);
m_theory = linspace(m_min, m_max, 1000);
B_theory = Gamma_half ./ ((m_theory - M).^2 + Gamma_half^2);
plot(m_theory, B_theory, 'r-', 'LineWidth', 2);
title('Theoretical Breit-Wigner Distribution');
xlabel('Mass (MeV)');
ylabel('B(m)');
grid on;

%% Statistical comparison
fprintf('Statistical Summary:\n');
fprintf('===================\n');
fprintf('Method (i) - Mean: %.2f MeV, Std: %.2f MeV\n', ...
    mean(accepted_masses_weighted), std(accepted_masses_weighted));
fprintf('Method (ii) - Mean: %.2f MeV, Std: %.2f MeV\n', ...
    mean(accepted_masses_reject), std(accepted_masses_reject));
fprintf('Theoretical peak at: %.1f MeV\n', M);

%% Additional comparison plot
figure('Position', [150, 150, 800, 600]);

% Overlay histograms
subplot(2, 1, 1);
[n1, x1] = hist(accepted_masses_weighted, 30);
[n2, x2] = hist(accepted_masses_reject, 30);

bar(x1, n1/max(n1), 'r');
hold on;
bar(x2, n2/max(n2), 'b');

legend('Weight Factor Method', 'Accept/Reject Method', 'Location', 'best');
title('Comparison of Both Methods (Normalized)');
xlabel('Mass (MeV)');
ylabel('Normalized Frequency');
grid on;

% Q-Q plot for comparison
subplot(2, 1, 2);
sorted1 = sort(accepted_masses_weighted);
sorted2 = sort(accepted_masses_reject);
plot(sorted1, sorted2, 'ko', 'MarkerSize', 3);
hold on;
min_val = min([sorted1; sorted2]);
max_val = max([sorted1; sorted2]);
plot([min_val, max_val], [min_val, max_val], 'r--');
title('Q-Q Plot: Method (i) vs Method (ii)');
xlabel('Weight Factor Method Quantiles');
ylabel('Accept/Reject Method Quantiles');
grid on;

