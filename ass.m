%%
%Q1
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
