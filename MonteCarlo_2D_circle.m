%% Monte Carlo simulation for photons in a circle

clear all;
close all;
clc;

% Parameters
R = 1; % Circle radius
tau = 10; % Optical depth
N_photons = 1e6; % Number of photons
max_scatterings = 300; % Maximum number of scatterings per photon
n_no_scatter_escape = 0; % Photons escaping without scattering
n_scatter_escape = zeros(1, max_scatterings); % Photons escaping after 1, 2, ... scatterings

% Mean free path
l_path = R / tau;

% Initial values
n_escaped = 0;
n_trapped = 0;
escaped_positions = [];
trapped_positions = [];

% Loop
for i = 1:N_photons
    % Initial position and direction
    position = [0, 0]; 
    direction = [0, -1]; % Moving along the +z axis

    for scattering_count = 0:max_scatterings
        % Distance until next interaction
        step_length = l_path * log(rand);
        
        new_position = position + step_length * direction;
        
        % Check if the photon escapes the circle
        if norm(new_position) > R
            n_escaped = n_escaped + 1;
            escaped_positions = [escaped_positions; new_position];
            if scattering_count == 0
                n_no_scatter_escape = n_no_scatter_escape + 1; % Escaped without scattering
            elseif scattering_count <= max_scatterings
                n_scatter_escape(scattering_count) = n_scatter_escape(scattering_count) + 1;
            end
            break; % Photon escapes, end simulation for this photon
        end
        
        % Photon remains inside, update position and scatter
        position = new_position;

        % Random direction (isotropic scattering)
        theta = 2 * pi * rand; % Random angle
        direction = [cos(theta), sin(theta)];
        
        if scattering_count == max_scatterings
            n_trapped = n_trapped + 1;
            trapped_positions = [trapped_positions; position];
        end
    end
end

% Percentages
P_escaped = (n_escaped / N_photons) * 100;
P_trapped = (n_trapped / N_photons) * 100;
P_no_scatter_escape = (n_no_scatter_escape / N_photons) * 100;
P_scatter_escape = (n_scatter_escape / N_photons) * 100;

P_theor_esc = exp(-tau) * 100;

% Results
fprintf('Total photons: %d (P = 100%%)\n', N_photons);
fprintf('Escaped photons: %d (P = %.2f%%)\n', n_escaped, P_escaped);
fprintf('Trapped photons: %d (P = %.2f%%)\n', n_trapped, P_trapped);
fprintf('Photons escaped without scattering: %d (P = %.2f%%)\n', n_no_scatter_escape, P_no_scatter_escape);
fprintf('Photons escaped theoretically without scattering P = %.2f%%)\n', P_theor_esc);

% Percentages for photons escaped after specific scatterings
for s = 1:max_scatterings
    fprintf('Photons escaped after %d scatterings: %d (P = %.2f%%)\n', s, n_scatter_escape(s), P_scatter_escape(s));
    if n_scatter_escape(s) == 0
        break;
    end
end

%% Plot

figure;
hold on;
rectangle('Position', [-R, -R, 2*R, 2*R], 'Curvature', [1, 1], 'EdgeColor', 'k', 'LineWidth', 5);

% Plot escaped photons
if ~isempty(escaped_positions)
    scatter(escaped_positions(:, 1), escaped_positions(:, 2), 2, 'm', 'filled'); % Escaped photons
end

% Plot trapped photons
if ~isempty(trapped_positions)
    scatter(trapped_positions(:, 1), trapped_positions(:, 2), 2, 'g', 'filled'); % Trapped photons
end
%xlim([-50 50]);
%ylim([-40 40]);
xlabel('x'); ylabel('z');
axis equal;
title('Photon Trajectories');
legend('Escaped Photons', 'Trapped Photons');
hold off;

%% Histogram

scatter_counts = 0:max_scatterings;
percentages = [P_no_scatter_escape, P_scatter_escape];

figure;
bar(scatter_counts, percentages, 'FaceColor', 'm');
hold on;
x = 3:0.1:45;
y = 15 * exp(- x/7);
%plot(x, y, 'k', 'LineWidth', 2);
xlim([0 300]);
xlabel('Number of Scatterings');
ylabel('Percentage of Escaped Photons (%)');
title(['Histogram of Escaped Photons by Number of Scatterings for ô = ' num2str(tau) ' ']);
grid on;
 