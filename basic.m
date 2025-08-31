%% Integrate 0 - 1 e^x dx with 1000 random numbers 

clc
clear all;

a=0;
b=1;
exact = exp(1)-1

sum = 0;
N = 1000;
h = (b-a)/N;

for i=1:N
    x = a + rand*(b-a);
    %% sum = sum + f(x)*h
    sum = sum + exp(x)*h;
end
sum

%% Radio active decay
clear all;
clc
N = 5000;
alpha = 0.03;
dt = 1;

for i = 0:dt:300
    k = 0;
    for j=1:N;
        if(rand<alpha*dt);
            k = k + 1;
        end 
    end
    N = N-k;
    x(i+1) = i;
    y(i+1) = N;
end

plot(x,y,'.');hold on;
N=5000;
for i = 0:dt:300
    z(i+1) = N*exp(-alpha*(i+1));
end
plot(x,z);hold off;


%% Hit and miss

clear all;
radius = 1;
count = 0;
for i=1:10000
    x(i) = rand*2-1;
    y(i) = rand*2-1;
    r = sqrt(x(i)^2 + y(i)^2);
    if(r<= radius)
        count = count+1;
        x2(count) = x(i);
        y2(count) = y(i);
    end
end
find_pi = 4*count/i
plot(x,y,'.',x2,y2,'r.');


%% MC by Direct sampling 
%% integration  0-1 dx/(1+x^2)
clear all;
clc; 

exact = 0.7850;    % Reference value
L = 20;            % Number of steps

z = zeros(1, L);         % Preallocate for efficiency
favg = zeros(1, L);      
sd = zeros(1, L);        

for k = 1:L
    N = 1000 * k;        % Sample size increases with k
    sumf = 0;            % Initialize sum
    sumf2 = 0;           % Initialize sum of squares
    for i = 1:N
        x = rand;                    % Random value in [0, 1]
        fx = 1 / (1 + (x^2));        % Function evaluation
        sumf = sumf + fx;            % Update sum
        sumf2 = sumf2 + fx^2;        % Update sum of squares
    end
    z(k) = k;                                    % Step count
    favg(k) = sumf / N;                          % Mean
    sd(k) = sqrt(((sumf2 / N) - (favg(k)^2)) / N);% Standard deviation over sqrt(N)
end

errorbar(z, favg, sd, 'o'); hold on;            % Error bars
plot(z, exact * ones(1, L), '-b'); hold off;     % Plot reference value across L
