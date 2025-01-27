%% Monte Carlo simulation for photons in a 3D sphere

clear all;
close all;
clc;

% Parameters
R = 1; % Sphere radius
tau = 3; % Optical depth
N_photons = 1e6; % Number of photons
max_scatterings = 100; % Maximum number of scatterings per photon
%absorption_prob = 0.3; % 30% probability of absorption

% Initialize counters
n_no_scatter_escape = 0; % Photons escaping without scattering
n_scatter_escape = zeros(1, max_scatterings); % Photons escaping after 1, 2, ... scatterings
%n_absorbed = 0; % Count absorbed photons
n_trapped = 0; % Count trapped photons

% Mean free path
l_path = R / tau;

% Initial values
n_escaped = 0;
escaped_positions = [];
trapped_positions = [];

% Loop
for i = 1:N_photons
    % Initial position and direction
    position = [0, 0, 0]; 
    direction = [0, 0, -1]; % Moving along the +z axis

    for scattering_count = 0:max_scatterings
        % Distance until next interaction
        step_length = l_path * log(rand);

        new_position = position + step_length * direction;

        % Check if the photon escapes the sphere
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

%         % Absorption check during scattering
%         if rand < absorption_prob
%             n_absorbed = n_absorbed + 1;
%             break; % Photon is absorbed and lost
%         end

        % Photon remains inside, update position and scatter
        position = new_position;

        % Random direction (isotropic scattering)
        theta = 2 * pi * rand; % Random azimuthal angle
        phi = acos(2 * rand - 1); % Random polar angle
        direction = [sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi)];

        if scattering_count == max_scatterings
            n_trapped = n_trapped + 1;
            trapped_positions = [trapped_positions; position]; 
        end
    end
end

% Calculate percentages
P_escaped = n_escaped / N_photons;
%P_absorbed = n_absorbed / N_photons;
P_trapped = n_trapped / N_photons;
P_no_scatter_escape = n_no_scatter_escape / N_photons;
P_scatter_escape = n_scatter_escape / N_photons;

% Results
fprintf('Total photons: %d (P = 100%%)\n', N_photons);
fprintf('Escaped photons: %d (P = %.2f%%)\n', n_escaped, P_escaped * 100);
%fprintf('Absorbed photons: %d (P = %.2f%%)\n', n_absorbed, P_absorbed * 100);
fprintf('Trapped photons: %d (P = %.2f%%)\n', n_trapped, P_trapped * 100);
fprintf('Photons escaped without scattering: %d (P = %.2f%%)\n', n_no_scatter_escape, P_no_scatter_escape * 100);

% Percentages for photons escaped after specific scatterings
for s = 1:max_scatterings
    fprintf('Photons escaped after %d scatterings: %d (P = %.2f%%)\n', s, n_scatter_escape(s), P_scatter_escape(s) * 100);
    if n_scatter_escape(s) == 0
        break;
    end
end

%% Visualization of photon escape and trapping
figure;
hold on;

[x, y, z] = sphere(50);
surf(R * x, R * y, R * z, 'FaceColor', 'm');

% Plot escaped photon positions
if ~isempty(escaped_positions)
    scatter3(escaped_positions(:, 1), escaped_positions(:, 2), escaped_positions(:, 3), ...
       1, 'b', 'filled', 'DisplayName', 'Escaped Photons'); % Escaped photons
end

% Plot trapped photon positions
if ~isempty(trapped_positions)
    scatter3(trapped_positions(:, 1), trapped_positions(:, 2), trapped_positions(:, 3), ...
        10, 'r', 'filled', 'DisplayName', 'Trapped Photons'); % Trapped photons
end

xlabel('x');
ylabel('y');
zlabel('z');
title('Photon Trajectories in a Sphere');
legend({'Sphere Boundary', 'Escaped Photons', 'Trapped Photons'});
axis equal;
grid on;
view(3);
hold off;

%% Histogram of escape percentages by scattering count
figure;
scatter_counts = 0:max_scatterings;
bar(scatter_counts, [P_no_scatter_escape, P_scatter_escape(1:max_scatterings)] * 100, 'FaceColor', [0, 0.447, 0.741]);
hold on;
x = 0:0.1:20;
y = 100 * exp(-1.8*x);
%plot(x, y, 'r', 'LineWidth', 2);
xlim([0 40]);
xlabel('Number of Scatterings');
ylabel('Percentage of Escaped Photons (%)');
title('Histogram of Escaped Photons by Scatterings');
grid on;
