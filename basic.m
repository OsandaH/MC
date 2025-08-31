%% Integrate 
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
    sum = sum + exp(x)*h;
end
sum
