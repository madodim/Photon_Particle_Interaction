clear all;
close all;
clc;

% Generate 1000 random numbers between 0 and 1
N = 100000;
random_numbers = rand(1, N);

% Statistics
mean_value = mean(random_numbers);
std_dev = std(random_numbers);

% Display statistics
fprintf('Generated %d random numbers.\n', N);
fprintf('Mean: %.4f\n', mean_value);
fprintf('Standard Deviation: %.4f\n', std_dev);

figure;
histogram(random_numbers, 10, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('Random Number Value');
ylabel('Frequency');
title('Histogram of Random Numbers');
grid on;

% Scatter plot
figure;
plot(random_numbers, '.', 'MarkerSize', 10);
xlabel('Index');
ylabel('Random Number Value');
title('Scatter Plot of Random Numbers');
grid on;
