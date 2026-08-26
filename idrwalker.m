function idrwalker(pdb_file, num_trials)
    % Automatically detect missing residue regions in a PDB file,
    % perform protein backbone Monte Carlo random walks, and check plausibility.
    % Based on the principles of the IDRWalker paper.

    %% 1. Read PDB and extract backbone atoms (N, CA, C)
    try
        pdb_struct = pdbread(pdb_file);
    catch ME
        fprintf('Error reading PDB file: %s\n', ME.message);
        return;
    end

    atoms = pdb_struct.Model.Atom;

    % Extract backbone atoms and their details
    backbone_atoms = struct('AtomName', {}, 'resSeq', {}, 'chainID', {}, 'X', {}, 'Y', {}, 'Z', {}); % Initialize as empty struct array with desired fields

    for i = 1:length(atoms)
        atom = atoms(i);
        % Check if the atom is a backbone atom (N, CA, or C)
        if strcmp(atom.AtomName, 'N') || strcmp(atom.AtomName, 'CA') || strcmp(atom.AtomName, 'C')
            % Create a new scalar structure with only the required fields
            new_backbone_atom = struct();
            new_backbone_atom.AtomName = atom.AtomName;
            new_backbone_atom.resSeq = atom.resSeq;
            new_backbone_atom.chainID = atom.chainID;
            new_backbone_atom.X = atom.X;
            new_backbone_atom.Y = atom.Y;
            new_backbone_atom.Z = atom.Z;

            % Append the new structure to the backbone_atoms array
            backbone_atoms = [backbone_atoms, new_backbone_atom];
        end
    end

    if isempty(backbone_atoms)
        fprintf('No backbone atoms (N, CA, C) found in the PDB file.\n');
        return;
    end

    %% 2. Identify missing residue regions (gaps in backbone)
    % Sort atoms by chain and residue number
    % Create a numerical or character representation for sorting chains
    chain_codes = cellfun(@double, {backbone_atoms.chainID}); % Convert chain IDs to numerical codes
    res_seqs = [backbone_atoms.resSeq];

    [~, sorted_indices] = sortrows([chain_codes', res_seqs']);
    backbone_atoms_sorted = backbone_atoms(sorted_indices);

    % Find unique residues and chains by creating formatted strings
    % Need to handle chainID as char

    % Create a cell array of strings, each combining ChainID and resSeq
    chain_res_strings = cell(length(backbone_atoms_sorted), 1);
    for k = 1:length(backbone_atoms_sorted)
        % Ensure chainID is a character array (it should be from PDBread, but good practice)
        chain_id_char = char(backbone_atoms_sorted(k).chainID);
        chain_res_strings{k} = sprintf('%s %d', chain_id_char, backbone_atoms_sorted(k).resSeq);
    end

    % Find unique combinations of chain and residue number
    unique_residues_cell = unique(chain_res_strings, 'stable'); % Use 'stable' to maintain original order

    % Initialize missing_regions BEFORE the loop that populates it
    missing_regions = {}; % Store information about missing regions

    for i = 1:size(unique_residues_cell, 1) - 1
        current_res_info = strtrim(unique_residues_cell{i,:});
        next_res_info = strtrim(unique_residues_cell{i+1,:});

        % Parse chainID and resSeq from the combined string
        % Find the space delimiter
        space_idx_current = find(isspace(current_res_info), 1, 'first');
        current_chain = current_res_info(1:space_idx_current-1);
        current_resnum = str2double(current_res_info(space_idx_current+1:end));

        space_idx_next = find(isspace(next_res_info), 1, 'first');
        next_chain = next_res_info(1:space_idx_next-1);
        next_resnum = str2double(next_res_info(space_idx_next+1:end));


        % Check for gap in residue number within the same chain
        if strcmp(current_chain, next_chain) && next_resnum > current_resnum + 1
            gap_start_res = current_resnum;
            gap_end_res = next_resnum;
            gap_len = gap_end_res - gap_start_res - 1;

            % Find the C, N, CA of the residue before the gap
            start_anchor_c = []; start_anchor_n = []; start_anchor_ca = [];
            end_anchor_c = []; end_anchor_n = []; end_anchor_ca = [];

            for j = 1:length(backbone_atoms_sorted)
                atom = backbone_atoms_sorted(j);
                if atom.resSeq == gap_start_res && strcmp(atom.chainID, current_chain)
                    if strcmp(atom.AtomName, 'C')
                        start_anchor_c = [atom.X, atom.Y, atom.Z];
                    elseif strcmp(atom.AtomName, 'N')
                        start_anchor_n = [atom.X, atom.Y, atom.Z];
                    elseif strcmp(atom.AtomName, 'CA')
                        start_anchor_ca = [atom.X, atom.Y, atom.Z];
                    end
                elseif atom.resSeq == gap_end_res && strcmp(atom.chainID, current_chain)
                     if strcmp(atom.AtomName, 'C')
                        end_anchor_c = [atom.X, atom.Y, atom.Z];
                    elseif strcmp(atom.AtomName, 'N')
                        end_anchor_n = [atom.X, atom.Y, atom.Z];
                    elseif strcmp(atom.AtomName, 'CA')
                        end_anchor_ca = [atom.X, atom.Y, atom.Z];
                    end
                end
            end

            % Ensure we found all necessary anchor atoms
            if ~isempty(start_anchor_c) && ~isempty(start_anchor_n) && ~isempty(start_anchor_ca) && ...
               ~isempty(end_anchor_c) && ~isempty(end_anchor_n) && ~isempty(end_anchor_ca)
                missing_regions{end+1} = struct(...
                    'chainID', current_chain, ...
                    'start_resnum', gap_start_res, ...
                    'end_resnum', gap_end_res, ...
                    'gap_len', gap_len, ...
                    'start_anchor_c', start_anchor_c, ...
                    'start_anchor_n', start_anchor_n, ...
                    'start_anchor_ca', start_anchor_ca, ...
                    'end_anchor_c', end_anchor_c, ...
                    'end_anchor_n', end_anchor_n, ...
                    'end_anchor_ca', end_anchor_ca);
            else
                 fprintf('Warning: Incomplete anchor atoms found for gap between residue %d (Chain %s) and %d. Skipping.\n', ...
                         gap_start_res, current_chain, gap_end_res);
            end
        end
    end

    if isempty(missing_regions)
        fprintf('No missing residue ranges detected based on backbone atoms.\n');
        return;
    end


    %% 3. For each missing region, run Monte Carlo random walks
    fprintf('Detected %d missing region(s) with complete backbone anchors:\n', length(missing_regions));

    % Subplot grid size
    ncols = 3; % Adjust as needed for desired plot layout
    nrows = ceil(length(missing_regions) / ncols);
    figure; % Create a single figure for all histograms

    for i = 1:length(missing_regions)
        region = missing_regions{i};
        fprintf('Gap %d: Chain %s, Residues %d–%d missing (%d residues)\n', ...
                i, region.chainID, region.start_resnum + 1, region.end_resnum - 1, region.gap_len);

        gap_len = region.gap_len;
        start_anchor_c = region.start_anchor_c;
        start_anchor_n = region.start_anchor_n;
        start_anchor_ca = region.start_anchor_ca;
        end_anchor_n = region.end_anchor_n; % We need the N of the residue AFTER the gap to connect to the last C

        % Target distance to connect the last generated C to the next N
        target_distance = norm(end_anchor_n - start_anchor_c);


        end_to_end_dist = zeros(num_trials, 1);

        % --- Monte Carlo Random Walks ---
        for t = 1:num_trials
            % Start the walk from the last C, N, CA of the preceding structured region
            prev_c_gen = start_anchor_c;
            prev_n_gen = start_anchor_n;
            prev_ca_gen = start_anchor_ca;

            simulation_successful = true;
            current_chain_c_coords = {prev_c_gen}; % Store generated C coordinates for this walk


            % Generate gap_len residues
            for j = 1:gap_len
                % Randomly sample dihedral angles (in radians)
                phi = -pi + (pi - (-pi)) * rand(); % Random float between -pi and pi
                psi = -pi + (pi - (-pi)) * rand(); % Random float between -pi and pi
                omega = pi; % Trans peptide bond

                % Generate coordinates for the next residue's backbone atoms
                [next_c, next_n, next_ca] = generateBackboneCoords(prev_c_gen, prev_n_gen, prev_ca_gen, phi, psi, omega);

                if isempty(next_c) % Check for failure (e.g., linear geometry)
                    simulation_successful = false;
                    break; % Stop generating residues for this walk
                end

                % Store the generated C coordinate
                current_chain_c_coords{j+1} = next_c;

                % Update previous atoms for the next iteration
                prev_c_gen = next_c;
                prev_n_gen = next_n;
                prev_ca_gen = next_ca;
            end

            % Calculate end-to-end distance if the walk was successful and has the correct length
            if simulation_successful && length(current_chain_c_coords) == gap_len + 1 % +1 for the starting C
                 % Distance between the last generated C and the target N anchor
                 last_generated_c = current_chain_c_coords{end};
                 distance = norm(end_anchor_n - last_generated_c);
                 end_to_end_dist(t) = distance;
            else
                 end_to_end_dist(t) = NaN; % Mark as failed simulation
            end
        end

        % Remove failed simulations
        end_to_end_dist = end_to_end_dist(~isnan(end_to_end_dist));

        % --- Analyze Results ---
        if ~isempty(end_to_end_dist)
            mu = mean(end_to_end_dist);
            sigma = std(end_to_end_dist);

            % Basic plausibility check: probability of being within a certain
            % range of the target distance under a normal distribution.
            % This is a simplification.
            % Check if sigma is effectively zero to avoid division by zero or Inf
            if sigma < 1e-9
                 if abs(target_distance - mu) < 0.5 % If target is very close to the single outcome
                     prob = 1;
                 else
                     prob = 0;
                 end
            else
                prob = normcdf(target_distance + 0.5, mu, sigma) - normcdf(target_distance - 0.5, mu, sigma);
            end

            % Report
            fprintf('  Target C-N Anchor distance = %.2f Å | Predicted μ = %.2f σ = %.2f | Plausibility P ≈ %.3f\n', ...
                target_distance, mu, sigma, prob);

            % Plot histogram
            subplot(nrows, ncols, i);
            histogram(end_to_end_dist, 40, 'FaceColor', [0.2 0.6 0.8]);
            hold on;
            ylims = ylim;
            plot([target_distance target_distance], ylims, 'r--', 'LineWidth', 1.5);
            title(sprintf('Gap %s:%d–%d (%d res)', region.chainID, region.start_resnum+1, region.end_resnum-1, gap_len), 'FontSize', 8);
            xlabel('End-to-End C-N Dist (Å)', 'FontSize', 7);
            ylabel('Frequency', 'FontSize', 7);
            set(gca, 'FontSize', 6);
            hold off;
        else
             fprintf('  No successful random walks for this gap.\n');
        end
    end
end

%% --- Helper Functions  ---

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

%% --- Helper Function: Rotate Vector ---
function rotated_vector = rotateVector(vector, axis, theta)
% rotateVector Rotates a vector around a given axis by an angle theta.
% Uses the rotationMatrix function.
    rot_mat = rotationMatrix(axis, theta);
    rotated_vector = (rot_mat * vector')'; % Apply rotation and transpose back
end

%% --- Helper Function: Rotation Matrix ---
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

