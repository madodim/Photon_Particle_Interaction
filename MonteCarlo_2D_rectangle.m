%% Monte Carlo simulation for photons in a rectangle

clear all;
close all;
clc;

% Parameters
D = 1; % Rectangle height (z-axis)
tau = 1; % Optical depth
N_photons = 1e6; % Number of photons
max_scatterings = 100; % Maximum number of scatterings per photon
n_no_scatter_escape = 0; % Photons escaping without scattering
n_scatter_escape = zeros(1, max_scatterings); % Photons escaping after 1, 2, ... scatterings

% Mean free path
l_path = D / tau;

% Initial values
n_escaped = 0;
n_trapped = 0;
n_reflected = 0; % Count of photons reflected (leaving from the bottom surface)
escaped_positions = [];
trapped_positions = [];
reflected_positions = [];

% Loop
for i = 1:N_photons
    % Initial position and direction
    position = [0, 0]; 
    direction = [0, -1]; % Moving along the +z axis

    for scattering_count = 0:max_scatterings
        % Distance until next interaction
        step_length = l_path * log(rand);
      
        new_position = position + step_length * direction;
        
        % Check if the photon escapes the rectangle through the upper side
        if new_position(2) >= D
            n_escaped = n_escaped + 1;
            escaped_positions = [escaped_positions; new_position];
            if scattering_count == 0
                n_no_scatter_escape = n_no_scatter_escape + 1; % Escaped without scattering
            elseif scattering_count <= max_scatterings
                n_scatter_escape(scattering_count) = n_scatter_escape(scattering_count) + 1;
            end
            break; % Photon escapes, end simulation for this photon
        end
        
        if new_position(2) <= 0
            n_reflected = n_reflected + 1;
            reflected_positions = [reflected_positions; new_position];
            break; % Photon reflects and exits, end simulation for this photon
        end

        % Photon remains inside, update position and scatter
        position = new_position;

        % Random direction (isotropic scattering in 2D)
        theta = 2 * pi * rand; % Random angle
        direction = [cos(theta), sin(theta)];
                
        % Check if the photon remains trapped
        if scattering_count == max_scatterings
            n_trapped = n_trapped + 1;
            trapped_positions = [trapped_positions; position];
        end    
    end
end

% Calculate percentages
P_escaped = n_escaped / N_photons;
P_trapped = n_trapped / N_photons;
P_reflected = n_reflected / N_photons;
P_no_scatter_escape = n_no_scatter_escape / N_photons;
P_scatter_escape = n_scatter_escape / N_photons; 

P_theor_esc = exp(-tau) * 100;

% Results
fprintf('Total photons: %d (P = 100%%)\n', N_photons);
fprintf('Escaped photons: %d (P = %.2f%%)\n', n_escaped, P_escaped * 100);
fprintf('Trapped photons: %d (P = %.2f%%)\n', n_trapped, P_trapped * 100);
fprintf('Reflected photons: %d (P = %.2f%%)\n', n_reflected, P_reflected * 100);
fprintf('Photons escaped without scattering: %d (P = %.2f%%)\n', n_no_scatter_escape, P_no_scatter_escape * 100);
fprintf('Photons escaped theoretically without scattering P = %.2f%%)\n', P_theor_esc);

% Percentages for photons escaped after specific scatterings
for s = 1:max_scatterings
    fprintf('Photons escaped after %d scatterings: %d (P = %.2f%%)\n', s, n_scatter_escape(s), P_scatter_escape(s) * 100);
    if n_scatter_escape(s) == 0
        break;
    end
end

%% Visualization of photon trajectories
figure;
hold on;

W = 200 * D;
rectangle('Position', [-W / 2, 0, W, D], 'EdgeColor', 'k', 'LineWidth', 2);

% Plot escaped photons
if ~isempty(escaped_positions)
    scatter(escaped_positions(:, 1), escaped_positions(:, 2), 2, 'm', 'filled'); % Escaped photons
end

% Plot trapped photons
if ~isempty(reflected_positions)
    scatter(reflected_positions(:, 1), reflected_positions(:, 2), 2, 'b', 'filled'); % Reflected photons
end

% Plot trapped photons
if ~isempty(trapped_positions)
    scatter(trapped_positions(:, 1), trapped_positions(:, 2), 2, 'b', 'filled'); % Trapped photons
end

xlabel('x');
ylabel('z');
xlim([-30 30]);
ylim([-20 20]);
title('Photon Trajectories in 2D Rectangle');
legend('Escaped Photons', 'Reflected Photons');
%axis equal;
grid on;
hold off;

%% Histogram of escaped photons by scattering number
scatter_counts = 0:max_scatterings; % Scatter numbers (0 to max_scatterings)
percentages = [P_no_scatter_escape, P_scatter_escape]; % Combine 0 scatter and others

figure;
bar(scatter_counts, percentages * 100, 'FaceColor', 'm');
hold on;
x = 0:0.1:20;
y = 25 * exp(- x/1.8);
plot(x, y, 'k', 'LineWidth', 2);

xlim([1 15]);
xlabel('Number of Scatterings');
ylabel('Percentage of Escaped Photons (%)');
title('Histogram of Escaped Photons by Number of Scatterings');
grid on;