%%%%%%%%%%%%% Q1_i
clc
clear all;

% Region: (x^2 + y^2 - 2x)^2 <= 4(x^2 + 4y^2)

N_val = [100,1000,10000];

xmin = -10;
xmax = 10;
ymin = -10;
ymax = 10;
box_area = (xmax-xmin)*(ymax-ymin);

for N=N_val
    x = xmin + (xmax-xmin)*rand(N,1)
    y = ymin + (ymax-ymin)*rand(N,1);

    inside = ((x.^2 + y.^2 -2*x).^2 <= 4*(x.^2 + 4*y.^2));

    area = box_area * sum(inside)/N;

    % Standard error
    p = sum(inside)/N;
    error = box_area * sqrt(p*(1-p)/N);

    fprintf('N=%d, Estimated Area = %.4f, Error ≈ %.4f\n', N, area, error);
    
    % Plot region with points
    figure;
    scatter(x(inside), y(inside)); hold on;
    scatter(x(~inside), y(~inside));
    title('Region (i): Monte Carlo Points');
    xlabel('x'); 
    ylabel('y'); 
    axis equal;
end

%%%%%%%%%%%%% Q1_ii
% Region: %
% -5/4 - 3/(3+x^2) <= y <= 5/4 + 3/(3+x^2)
% -1/3 - 5*y^2/6 + y^4/6 <= x <= 11/3 - (2/3)*sqrt(2)*|y|^(3/2)

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

%%%%%%%%%%%%% Q2_MC

clc
clear all;

f = @(x) exp(-x.^2);

exact = integral(f,0,1);
disp("Exact Value - "+ exact);

N = 1000;
x_rand = rand(N,1);
I = mean(f(x_rand));

error = abs(I-exact);

fprintf('Plain Monte Carlo:\n');
fprintf('Estimate = %.8f, Error = %.2e\n', I, error);

%%%%%%%%%%%%% Q2_Stratified_Sampling
clc; 
clear; 
close all;

f = @(x) exp(-x.^2);
true_val = integral(f, 0, 1)

N = 1000; 
edges = linspace(0,1,N+1);
strat_samples = zeros(N,1);

for i = 1:N
    u = rand;
    x_strat = edges(i) + u*(edges(i+1)-edges(i));
    strat_samples(i) = f(x_strat);
end
I = mean(strat_samples);
error = abs(I-true_val);
fprintf('Stratified Sampling:\n');
fprintf('Estimate = %.8f, Error = %.2e\n', I, error);

