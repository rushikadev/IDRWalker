% Script to study the protein backbone random walk law using Monte Carlo simulation
% This script simulates random walks of different lengths to analyze
% the distribution of end-to-end distances and their relationship with
% the number of residues, as described in the IDRWalker paper's study
% of the random walk law.
% Helper functions generateBackboneCoords, rotateVector, and rotationMatrix
% are included as local functions within this file.

% Simulation parameters
num_simulations_per_length = 1000; % Number of Monte Carlo trials for each chain length

% Define different chain lengths (number of residues) to simulate
% Include small and large numbers as described in the paper
chain_lengths = [5, 10, 20, 50, 100, 150, 200]; % Example lengths

% Initialize storage for results
mean_distances = zeros(size(chain_lengths));
std_distances = zeros(size(chain_lengths));
all_end_to_end_distances = cell(size(chain_lengths)); % To store all distances for histograms

fprintf('Starting Monte Carlo simulations for different chain lengths...\n');

% --- Run Simulations for Each Chain Length ---
for i = 1:length(chain_lengths)
    num_residues_in_walk = chain_lengths(i);
    fprintf('Simulating walks of %d residues (%d trials)...\n', num_residues_in_walk, num_simulations_per_length);

    end_to_end_distances_for_length = zeros(num_simulations_per_length, 1);
    successful_simulations = 0;

    % --- Monte Carlo Simulation Loop (Core Random Walk Logic) ---
    % Initialize the first few atoms to start the chain (idealized geometry)
    % N0 at (0,0,0), Ca0 along x-axis, and C0 in xy-plane
    l_n_ca = 1.458;
    l_ca_c = 1.525;
    angle_n_ca_c = deg2rad(111.2); % N-Ca-C angle

    first_n = [0.0, 0.0, 0.0];
    first_ca = [l_n_ca, 0.0, 0.0];
    x_c0 = first_ca(1) + l_ca_c * cos(pi - angle_n_ca_c);
    y_c0 = first_ca(2) + l_ca_c * sin(pi - angle_n_ca_c);
    first_c = [x_c0, y_c0, 0.0];

    for t = 1:num_simulations_per_length
        prev_c_gen = first_c;
        prev_n_gen = first_n;
        prev_ca_gen = first_ca;

        simulation_successful = true;

        % We generate num_residues_in_walk - 1 additional residues
        for j = 1:(num_residues_in_walk - 1)
            % Randomly sample dihedral angles (in radians)
            % Phi and Psi ranges are typically -pi to +pi
            % Omega is usually pi radians (trans peptide bond)
            phi = -pi + (pi - (-pi)) * rand();
            psi = -pi + (pi - (-pi)) * rand();
            omega = pi; % Trans peptide bond

            % Generate coordinates for the next residue's backbone atoms using local function
            [next_c, next_n, next_ca] = generateBackboneCoords(prev_c_gen, prev_n_gen, prev_ca_gen, phi, psi, omega);

            if isempty(next_c) % Check for failure (e.g., linear geometry)
                 simulation_successful = false;
                 break; % Stop this simulation
            end

            % Update previous atoms for the next iteration
            prev_c_gen = next_c;
            prev_n_gen = next_n;
            prev_ca_gen = next_ca;
        end

        % Calculate end-to-end distance if the full chain was generated successfully
        if simulation_successful
            % End-to-end distance: distance between the first N and the last C
            start_point = first_n; % Starting point is the first N atom
            end_point = prev_c_gen; % Ending point is the last generated C atom
            distance = norm(end_point - start_point);
            successful_simulations = successful_simulations + 1;
            end_to_end_distances_for_length(successful_simulations) = distance;
        end
    end

    % Store results for this chain length
    end_to_end_distances_for_length = end_to_end_distances_for_length(1:successful_simulations);
    if ~isempty(end_to_end_distances_for_length)
        mean_distances(i) = mean(end_to_end_distances_for_length);
        std_distances(i) = std(end_to_end_distances_for_length);
        all_end_to_end_distances{i} = end_to_end_distances_for_length;
    else
        mean_distances(i) = NaN;
        std_distances(i) = NaN;
        all_end_to_end_distances{i} = [];
    end

    fprintf('  Completed %d successful simulations.\n', successful_simulations);
end

fprintf('Simulations finished.\n');

% --- Analysis and Visualization ---

% 1. Show distribution for small and large number of steps
figure;

% Find indices for a small and a large number of steps
% You can adjust these indices to pick different lengths from chain_lengths
[~, small_idx] = min(abs(chain_lengths - 10)); % Example: find length closest to 10
[~, large_idx] = min(abs(chain_lengths - 100)); % Example: find length closest to 100

% Plot histogram for a small number of steps
subplot(1, 2, 1);
histogram(all_end_to_end_distances{small_idx}, 'FaceColor', [0.2 0.6 0.8]);
title(sprintf('End-to-End Distance Distribution (%d residues)', chain_lengths(small_idx)));
xlabel('Distance (Å)');
ylabel('Frequency');
grid on;

% Plot histogram and normal distribution overlay for a large number of steps
subplot(1, 2, 2);
distances_large = all_end_to_end_distances{large_idx};
if ~isempty(distances_large)
    histogram(distances_large, 'FaceColor', [0.2 0.6 0.8], 'Normalization', 'pdf'); % Use pdf normalization for comparison with normal distribution
    hold on;

    % Fit a normal distribution to the data
    pd = fitdist(distances_large, 'Normal');

    % Plot the fitted normal distribution
    x_values = linspace(min(distances_large), max(distances_large), 100);
    plot(x_values, pdf(pd, x_values), 'r-', 'LineWidth', 2);
    hold off;

    title(sprintf('End-to-End Distance Distribution (%d residues) with Normal Fit', chain_lengths(large_idx)));
    xlabel('Distance (Å)');
    ylabel('Probability Density');
    legend('Simulation Data', 'Normal Distribution Fit');
    grid on;
else
    title(sprintf('No successful simulations for %d residues to plot distribution', chain_lengths(large_idx)));
end


% 2. Show relationship between number of steps and mean distance
figure;
plot(chain_lengths, mean_distances, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;

% Perform linear fitting (mean distance vs. number of steps) - as described in the paper
% Note: For an ideal random coil, the RMS end-to-end distance is proportional to sqrt(N)
% However, the paper mentions linear fitting of "the relationship between ... and the number of steps".
% Let's fit mean distance vs. number of steps as a simple interpretation.
% A more theoretically correct fit would be RMS distance vs. sqrt(N).

% Remove NaN values from fitting if some simulations failed
valid_lengths = chain_lengths(~isnan(mean_distances));
valid_mean_distances = mean_distances(~isnan(mean_distances));

if length(valid_lengths) >= 2
    p = polyfit(valid_lengths, valid_mean_distances, 1); % Fit a 1st degree polynomial (linear)
    y_fit = polyval(p, valid_lengths);
    plot(valid_lengths, y_fit, 'r--', 'LineWidth', 1.5);
    legend('Simulated Mean Distance', sprintf('Linear Fit (y = %.2fx + %.2f)', p(1), p(2)));

    fprintf('\nLinear fit results (Mean Distance vs. Number of Steps):\n');
    fprintf('  Slope: %.4f\n', p(1));
    fprintf('  Intercept: %.4f\n', p(2));

    % Optional: Fit RMS distance vs. sqrt(N) (more theoretically accurate for random coils)
    % Calculate RMS distance for each length from all successful simulations
    rms_distances = zeros(size(chain_lengths));
    for k = 1:length(chain_lengths)
        if ~isempty(all_end_to_end_distances{k})
            rms_distances(k) = sqrt(mean(all_end_to_end_distances{k}.^2));
        else
            rms_distances(k) = NaN;
        end
    end

    valid_rms_distances = rms_distances(~isnan(rms_distances));
    sqrt_valid_lengths = sqrt(chain_lengths(~isnan(rms_distances)));


    if length(sqrt_valid_lengths) >= 2
        p_rms = polyfit(sqrt_valid_lengths, valid_rms_distances, 1);
        y_fit_rms = polyval(p_rms, sqrt_valid_lengths);
        figure; % New figure for RMS vs sqrt(N)
        plot(sqrt_valid_lengths, valid_rms_distances, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
        hold on;
        plot(sqrt_valid_lengths, y_fit_rms, 'g--', 'LineWidth', 1.5);
        legend('Simulated RMS Distance', sprintf('Linear Fit (y = %.2fx + %.2f)', p_rms(1), p_rms(2)));
        xlabel('Square Root of Number of Residues (\surd{N})');
        ylabel('RMS End-to-End Distance (Å)');
        title('Relationship between RMS Distance and \surd{N}');
        grid on;

        fprintf('\nLinear fit results (RMS Distance vs. sqrt(N)):\n');
        fprintf('  Slope (Persistence Length estimate): %.4f\n', p_rms(1));
        fprintf('  Intercept: %.4f\n', p_rms(2));
        fprintf('Note: For an ideal random coil, the RMS end-to-end distance is proportional to the square root of the number of steps.\n');
    end

else
    fprintf('\nNot enough data points to perform linear fitting.\n');
end


xlabel('Number of Residues');
ylabel('Mean End-to-End Distance (Å)');
title('Relationship between Mean End-to-End Distance and Number of Residues');
grid on;
hold off;


%% --- Local Helper Functions ---

function [next_c, next_n, next_ca] = generateBackboneCoords(prev_c, prev_n, prev_ca, phi, psi, omega)
% generateBackboneCoords Generates coordinates for the next set of backbone
% atoms using idealized geometry and dihedral angles.

% Idealized bond lengths (Angstroms)
l_n_ca = 1.458;
l_ca_c = 1.525;
l_c_n = 1.329;

% Idealized bond angles (radians)
angle_c_n_ca = deg2rad(121.7); % C-(next N)-Ca
angle_n_ca_c = deg2rad(111.2); % N-Ca-C
angle_ca_c_n = deg2rad(116.2); % Ca-C-(next N)

% Step 1: Generate next N from prev_c
% r1 = prev_c. We are generating r2 = next_n.
% v1 is vector from prev_ca to prev_c
% n1 is normal to plane defined by prev_n, prev_ca, prev_c

vector1_plane = prev_ca - prev_n;
vector2_plane = prev_c - prev_ca;
n1_for_next_n = cross(vector1_plane, vector2_plane);
if norm(n1_for_next_n) < 1e-9 % Use a small tolerance for floating point comparison
     % Handle linear case
     next_c = []; next_n = []; next_ca = []; % Indicate failure
     return;
end
n1_for_next_n = n1_for_next_n / norm(n1_for_next_n);

v1_ca_c = prev_c - prev_ca;
if norm(v1_ca_c) < 1e-9 % Handle case where prev_ca and prev_c are the same
     next_c = []; next_n = []; next_ca = []; % Indicate failure
     return;
end
v1_ca_c_norm = v1_ca_c / norm(v1_ca_c);

% theta_d is omega (angle around Ca-C bond)
theta_d_omega = omega;
% theta_a is angle_ca_c_n (angle at C)
theta_a_ca_c_n_rad = angle_ca_c_n;

n2_for_next_n = rotateVector(n1_for_next_n, v1_ca_c_norm, theta_d_omega);
v2_c_nextn = rotateVector(-v1_ca_c_norm, n2_for_next_n, theta_a_ca_c_n_rad);

next_n = prev_c + l_c_n * v2_c_nextn;


% Step 2: Generate next Ca from next_n
% r1 = next_n. We are generating r2 = next_ca.
% v1 is along C-N
% n1 is normal to the plane defined by prev_ca, prev_c, next_n.

vector1_plane = prev_c - prev_ca;
vector2_plane = next_n - prev_c;
n1_for_next_ca = cross(vector1_plane, vector2_plane);
if norm(n1_for_next_ca) < 1e-9
     next_c = []; next_n = []; next_ca = []; % Indicate failure
     return;
end
n1_for_next_ca = n1_for_next_ca / norm(n1_for_next_ca);

v1_c_nextn = next_n - prev_c;
if norm(v1_c_nextn) < 1e-9
     next_c = []; next_n = []; next_ca = []; % Indicate failure
     return;
end
v1_c_nextn_norm = v1_c_nextn / norm(v1_c_nextn);

% theta_d is phi (angle around C-N bond)
theta_d_phi = phi;
% theta_a is angle_c_n_ca (angle at N)
theta_a_c_n_ca_rad = angle_c_n_ca;

n2_for_next_ca = rotateVector(n1_for_next_ca, v1_c_nextn_norm, theta_d_phi);
v2_nextn_nextca = rotateVector(-v1_c_nextn_norm, n2_for_next_ca, theta_a_c_n_ca_rad);

next_ca = next_n + l_n_ca * v2_nextn_nextca;

% Step 3: Generate next C from next_ca
% r1 = next_ca. We are generating r2 = next_c.
% v1 is along N-Ca
% n1 is normal to the plane defined by prev_c, next_n, next_ca.

vector1_plane = next_n - prev_c;
vector2_plane = next_ca - next_n;
n1_for_next_c = cross(vector1_plane, vector2_plane);
if norm(n1_for_next_c) < 1e-9
     next_c = []; next_n = []; next_ca = []; % Indicate failure
     return;
end
n1_for_next_c = n1_for_next_c / norm(n1_for_next_c);

v1_nextn_nextca = next_ca - next_n;
if norm(v1_nextn_nextca) < 1e-9
     next_c = []; next_n = []; next_ca = []; % Indicate failure
     return;
end
v1_nextn_nextca_norm = v1_nextn_nextca / norm(v1_nextn_nextca);

% theta_d is psi (angle around N-Ca bond)
theta_d_psi = psi;
% theta_a is angle_n_ca_c (angle at Ca)
theta_a_n_ca_c_rad = angle_n_ca_c;

n2_for_next_c = rotateVector(n1_for_next_c, v1_nextn_nextca_norm, theta_d_psi);
v2_nextca_nextc = rotateVector(-v1_nextn_nextca_norm, n2_for_next_c, theta_a_n_ca_c_rad);

next_c = next_ca + l_ca_c * v2_nextca_nextc;

end

%% --- Local Helper Function: Rotate Vector ---
function rotated_vector = rotateVector(vector, axis, theta)
% rotateVector Rotates a vector around a given axis by an angle theta.
% Uses the rotationMatrix function.
    rot_mat = rotationMatrix(axis, theta);
    rotated_vector = (rot_mat * vector')'; % Apply rotation and transpose back
end

%% --- Local Helper Function: Rotation Matrix ---
function R = rotationMatrix(axis, theta)
% rotationMatrix Calculates the rotation matrix associated with
% counterclockwise rotation about the given axis by theta radians.
    axis = axis / norm(axis);
    a = cos(theta / 2.0);
    b = -axis(1) * sin(theta / 2.0);
    c = -axis(2) * sin(theta / 2.0);
    d = -axis(3) * sin(theta / 2.0);
    R = [a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c);
         2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b);
         2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c];
end
