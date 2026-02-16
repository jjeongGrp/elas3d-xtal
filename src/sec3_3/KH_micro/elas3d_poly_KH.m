%% Synthetic 3D Microstructure Generator for ELAS3D-Xtal
% ------------------------------------------------------------------------------
% This script generates realistic 3D voxelized polycrystals containing
% columnar or equiaxed grains alongside process-induced defects (pores).
% Grains are seeded using statistical distributions (size and aspect ratio),
% and experimental defect data can be mapped directly onto the mesh.
% The final voxel phase map and crystallographic orientations are exported
% to an HDF5 file for direct ingestion into the ELAS3D-Xtal solver.
%
% Dependency:
%   - MTEX Toolbox (required for orientation assignment and IPF coloring)
%
% Inputs:
%   - Script configuration parameters (Domain size, Grid resolution, Statistics)
%   - Defect data file (Excel/CSV) containing pore centers (x,y,z) and diameters
%
% Outputs:
%   - 'input_structure_poly.h5': HDF5 file containing the 3D phase map and orientations
%   - Workspace data (.mat) and visualization figures (Cross-sections, 3D renders)
%
% Core Workflow:
%   1. Define computational domain and voxel discretization.
%   2. Seed grains based on lognormal size and aspect ratio distributions.
%   3. Assign grain phases via parallelized spatial filtering.
%   4. Map explicit defect/pore data into the voxel grid.
%   5. Assign crystallographic orientations (Random or Fiber texture).
%
% Authors:       Juyoung Jeong and Veera Sundararaghavan
% Affiliation:   Department of Aerospace Engineering, University of Michigan,
%                Ann Arbor, MI 48109, USA.
% Repository:    https://github.com/jjeongGrp/elas3d-xtal
%
% Acknowledgement:
%   This work was supported by the Defense Advanced Research Projects Agency
%   (DARPA) SURGE program under Cooperative Agreement No. HR0011-25-2-0009,
%   'Predictive Real-time Intelligence for Metallic Endurance (PRIME)'.
%
% License:       MIT License (See LICENSE file in the repository)
%
% ------------------------------------------------------------------------------
%% Columnar/Equiaxed Grain Microstructure with Lognormal Aspect Ratios
format compact; clear; close all;
tStart = tic;

%% 1. Configuration & Parameters -------------------------------------------------------------
% --- Domain Variables ---
max_len = 5.0;     
cube_len = 5.0;     
grid_num = 800;     
dx = cube_len/grid_num;
xmin = -2.5;
ymin = -2.5;
zmin = -2.5;

% --- Grain Statistics Variables ---
mean_d = 0.0565;       
std_d  = 0.25*mean_d; 

% Aspect ratio variables
mean_aspect = 3.253;  
std_aspect = 0.25 * mean_aspect;

% --- Orientation Settings Variables ---
orientation_type = 'textured'; 
sigma_spread = 12.5;         

% --- General Descriptor Variables ---
pore_type = 'KH';
desc = 'Keyhole Voids';

% --- Output Settings Variables ---
save_data = 'off';
elas3d_input = 'off';
phase_data_type = 'uint32';


% --- Visualization Variables -
stride = 1; 
microstructure3d_type = 'surface';
max_pores_show = 1;
num_grains_iso = 10;
micro_plot = 'on';
K_grain = 2000; 

% --- Pore Settings Variables ---
pore = 'on';
input_pore_file = 'Defect_Locations_KH.xlsx';

%% 2. Pore Data -----------------------------------------------------------

if strcmp(pore, 'on')
    % Load the entire defect dataset into a matrix
    data = readmatrix(input_pore_file);

    % Extract spatial coordinates of pore centers
    x_pore       = data(:,9);      % [mm] x_pore: X-coordinate of the pore center of mass.
    y_pore       = data(:,10);     % [mm] y_pore: Y-coordinate of the pore center of mass.
    z_pore       = data(:,11);     % [mm] z_pore: Z-coordinate of the pore center of mass.

    % Extract size metrics
    d = data(:,8);                 % [mm] d: Equivalent spherical diameter of the pore.
    r = d/2;                       % [mm] r: Calculated radius of the pore.
    V_pore = data(:,7);            % [mm³] V_pore: Physical volume of the pore.
    D_Feret_max = data(:,12);      % [mm] D_Feret_max: Maximum distance between any two points on the pore boundary.

    % Count total number of defects to inject
    N_pores = numel(x_pore);       % [integer] N_pores: Total number of pores loaded from the file.

elseif strcmp(pore, 'off')
    % Fallback if the specified excel file does not exist
    x_pore = [];
    y_pore = [];
    z_pore = [];
    d = [];
    r = [];
    V_pore = [];
    D_Feret_max = [];
    N_pores = 0;
    warning('No pore file found. Proceeding with no pores.');
end

%% 3. Computational Domain ------------------------------------------------

% Enforce max length constraint
cube_len = min(cube_len, max_len); % [mm] Safety check, though not strictly needed

% --- Define Domain Boundaries ---
xmax = xmin + cube_len; % [mm] xmax: Upper boundary of the domain on the X-axis.
ymax = ymin + cube_len; % [mm] ymax: Upper boundary of the domain on the Y-axis.
zmax = zmin + cube_len; % [mm] zmax: Upper boundary of the domain on the Z-axis.

% [array] cuboid_size: 1x3 vector containing the absolute lengths of the domain in X, Y, and Z.
cuboid_size = [xmax-xmin, ymax-ymin, zmax-zmin];

% --- Generate 1D and 3D Grids ---
xvec = xmin:dx:xmax;    % [array] xvec: 1D array of X-coordinates spaced by voxel size.
Nx = length(xvec);      % [integer] Nx: Number of grid nodes along the X-axis.

yvec = ymin:dx:ymax;    % [array] yvec: 1D array of Y-coordinates spaced by voxel size. 
Ny = length(yvec);      % [integer] Ny: Number of grid nodes along the Y-axis.

zvec = zmin:dx:zmax;    % [array] zvec: 1D array of Z-coordinates spaced by voxel size.
Nz = length(zvec);      % [integer] Nz: Number of grid nodes along the Z-axis.

% [3D matrices] X, Y, Z: Full 3D spatial coordinate arrays generated from 1D vectors.
[X, Y, Z] = ndgrid(xvec, yvec, zvec);


% [mm³] domain_vol: Total physical volume of the computational domain.
domain_vol = prod(cuboid_size);


%% 4. Lognormal Grain Size and Aspect Ratio --------------------------------

% Estimate grain count for target volume & mean diameter (cylinder approx)
% [mm³] mean_grain_vol: Theoretical average volume of a single grain, assuming a cylindrical shape (Area * Height).
mean_grain_vol = pi*(mean_d/2)^2 * (mean_aspect*mean_d);

% [integer] N_grains: Estimated total number of grains needed to pack the entire simulation domain.
N_grains = ceil(domain_vol / mean_grain_vol);

% --- Grain sizes (diameter in XY) ---
% Convert standard mean and std dev into lognormal distribution parameters
sigma_log = sqrt(log((std_d/mean_d)^2 + 1)); % [dimensionless] Lognormal shape parameter (standard deviation of the log of the distribution).

mu_log = log(mean_d) - 0.5*sigma_log^2;  % [dimensionless] Lognormal scale parameter (mean of the log of the distribution).

grain_diameters = lognrnd(mu_log, sigma_log, N_grains, 1); % [array] Nx1 array of randomly sampled diameters following the specified lognormal curve [mm].

r_xy = grain_diameters / 2;   % [array] Nx1 array representing the primary radius of each grain in the transverse (XY) plane [mm].

% --- Aspect ratio per grain (length/width) ---
% Convert standard mean and std dev of aspect ratios into lognormal parameters
sigma_log_a = sqrt(log((std_aspect/mean_aspect)^2 + 1)); % [dimensionless] Lognormal shape parameter for aspect ratios.

mu_log_a    = log(mean_aspect) - 0.5*sigma_log_a^2; % [dimensionless] Lognormal scale parameter for aspect ratios.

aspect_ratios = lognrnd(mu_log_a, sigma_log_a, N_grains, 1); % [array][dimensionless] Nx1 array of randomly sampled elongation ratios (Z-length / XY-width).

r_z  = aspect_ratios .* r_xy; % [array] Nx1 array representing the radius (half-length) of each grain along the vertical (Z) axis [mm].

% --- Grain seeds (random in the box) ---
% Generate random XYZ coordinate triplets within the bounding box of the domain
grain_seeds = [rand(N_grains,1)*cuboid_size(1)+xmin, ...
    rand(N_grains,1)*cuboid_size(2)+ymin, ...
    rand(N_grains,1)*cuboid_size(3)+zmin]; % [matrix] N_grains x 3 array containing the (X, Y, Z) center points for each grain seed [mm].

%% 5.Anisotropic Assignment -------------------------------------

% 5-1: Print header for upcoming block assignment step
fprintf('\n=== PARALLEL + SPATIAL FILTERING MICROSTRUCTURE GENERATION ===\n');

main_timer = tic; % [timer] Tracks the total time taken for the entire assignment phase.


% Setup parallel pool
if isempty(gcp('nocreate'))
    fprintf('Starting parallel pool...\n');
	
    pool_timer = tic; % [timer] Tracks how long it takes to spin up the parallel workers.
	
	% [integer] num_workers: Automatically detects the number of physical CPU cores available.
    num_workers = feature('numcores');
	
	% Initialize the parallel pool using local workers.
    parpool('local', num_workers);
    fprintf('Parallel pool started in %.1fs with %d workers\n', toc(pool_timer), gcp().NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', gcp().NumWorkers);
end

% Pre-allocate memory for phase assignment and computed properties
phase = zeros(Nx, Ny, Nz, phase_data_type); % [3D matrix] The primary output grid (Nx x Ny x Nz) that will store the integer Grain ID for every voxel.

% --- Memory Optimizations ---
% Pre-calculate the inverse of the squared radii to replace slow division with fast multiplication later.
inv_r_xy2 = single(1 ./ (r_xy.^2)); % [array] Nx1 array of 1/(r_xy^2) cast to single precision to save memory.

inv_r_z2 = single(1 ./ (r_z.^2)); % [array] Nx1 array of 1/(r_z^2) cast to single precision.

grain_seeds_single = single(grain_seeds); % [matrix] Copy of grain coordinates cast to single precision to reduce memory bandwidth in parallel loops.

% 5-2 Build Spatial Acceleration Structure
fprintf('\nBuilding spatial acceleration structure...\n');

accel_timer = tic; % [timer] Tracks the time taken to build the spatial search grid.

grid_res = 15; % [integer] Number of coarse subdivisions along each axis for the spatial filter.

% [array] x_edges, y_edges, z_edges: 1D vectors defining the boundaries of the coarse macro-cells.
x_edges = linspace(min(xvec), max(xvec), grid_res+1); 
y_edges = linspace(min(yvec), max(yvec), grid_res+1);
z_edges = linspace(min(zvec), max(zvec), grid_res+1);

% Create spatial grid using cell array indexing that parfor can handle
% Compute maximum possible influence region for any grain
max_influence = max([r_xy(:); r_z(:)]) * 1.5; % [float] The absolute maximum distance any grain can reach (uses the largest radius found, plus a 50% safety buffer).

total_cells = grid_res^3; % [integer] Total number of macro-cells in the acceleration grid.

% Prepare output acceleration grid as a cell array (linear indexing for parfor compatibility)
% Pre-allocate output cell array
spatial_grid_temp = cell(total_cells, 1); % [cell array] Stores a list of relevant Grain IDs for each of the macro-cells using 1D indexing.

% Parallel construction using linear indexing
% For each cell, find which grains may influence it
parfor linear_idx = 1:total_cells

    % Convert the 1D linear index back into 3D subscript indices (i, j, k) for the coarse grid.
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);

    % Define the physical bounding box [min, max] for the current macro-cell.
    cell_xmin = x_edges(i); cell_xmax = x_edges(i+1);
    cell_ymin = y_edges(j); cell_ymax = y_edges(j+1);
    cell_zmin = z_edges(k); cell_zmax = z_edges(k+1);

    % Find all grains whose centers are close enough to potentially affect this cell
    relevant_grains = find(...
        grain_seeds_single(:,1) >= (cell_xmin - max_influence) & ...
        grain_seeds_single(:,1) <= (cell_xmax + max_influence) & ...
        grain_seeds_single(:,2) >= (cell_ymin - max_influence) & ...
        grain_seeds_single(:,2) <= (cell_ymax + max_influence) & ...
        grain_seeds_single(:,3) >= (cell_zmin - max_influence) & ...
        grain_seeds_single(:,3) <= (cell_zmax + max_influence));     % [array] A filtered list of Grain IDs that fall within the expanded bounding box of this macro-cell.

    % Store the filtered list of relevant grains into the current cell of the acceleration structure.
    spatial_grid_temp{linear_idx} = relevant_grains;
end

% Convert back to 3D indexing for easier access
spatial_grid = cell(grid_res, grid_res, grid_res); % [cell array] Converts the 1D temp array into a 3D cell array for easier (i,j,k) lookups later.
for linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);
    spatial_grid{i,j,k} = spatial_grid_temp{linear_idx};
end
clear spatial_grid_temp; % Free memory

fprintf('Spatial acceleration structure built in %.1fs\n', toc(accel_timer));

% 5-3 Setup Parallel Block Processing

% Set block size for actual parallel assignment loop
block_size = 50; % [integer] Number of voxels along each edge of a processing block (e.g., a 50x50x50 voxel chunk).

% [integer] n_blocks_x/y/z: Calculates how many blocks fit along each axis.
n_blocks_x = ceil(Nx / block_size);
n_blocks_y = ceil(Ny / block_size);
n_blocks_z = ceil(Nz / block_size);

total_blocks = n_blocks_x * n_blocks_y * n_blocks_z; % [integer] Total number of discrete chunks the domain is divided into.

% Report block setup
fprintf('Processing %d blocks (%dx%dx%d) in parallel:\n', total_blocks, n_blocks_x, n_blocks_y, n_blocks_z);

% Create block processing batches
% Calculate how to group the blocks based on available CPU cores to maximize worker uptime.
num_workers = gcp().NumWorkers;

blocks_per_batch = max(1, floor(total_blocks / (num_workers * 4))); % [integer] Number of blocks grouped together to be sent to workers at once.

n_batches = ceil(total_blocks / blocks_per_batch); % [integer] Total number of iterations required to process all blocks.

fprintf('Using %d batches of %d blocks each for progress tracking\n', n_batches, blocks_per_batch);

% Initialize progress tracking variables
completed_blocks = 0;  % [integer] Running tally of blocks successfully assigned.

total_grains_processed = 0; % [integer] Running tally of how many individual grain distance calculations occurred.

blocks_with_grains = 0; % [integer] Tally of blocks that actually contained structural data (weren't empty space).

batch_times = []; % [array] Array storing the time taken for each batch to compute ETA.

% 5-4 Process Batches with Progress Tracking
fprintf('\nStarting parallel processing with real-time progress:\n');

process_timer = tic; % [timer] Tracks total execution time of the block processing phase.

for batch_idx = 1:n_batches
    batch_timer = tic; % [timer] Tracks execution time for the current batch.

    % Determine the start and end indices for blocks in the current batch
    batch_start = (batch_idx - 1) * blocks_per_batch + 1;
    batch_end = min(batch_idx * blocks_per_batch, total_blocks);
    current_batch_size = batch_end - batch_start + 1;

    % Create list of block indices for current batch
    block_list_batch = zeros(current_batch_size, 4); % [matrix] Stores the [X, Y, Z, Linear] block IDs assigned to this batch.
    for local_idx = 1:current_batch_size
        block_num = batch_start + local_idx - 1;
		% Convert linear block ID back to 3D coordinates (bz, by, bx)
        bz = ceil(block_num / (n_blocks_x * n_blocks_y));
        remaining = block_num - (bz-1) * n_blocks_x * n_blocks_y;
        by = ceil(remaining / n_blocks_x);
        bx = remaining - (by-1) * n_blocks_x;
        block_list_batch(local_idx, :) = [bx, by, bz, block_num];
    end

    % Pre-allocate cell arrays to collect results from parallel workers
    batch_results = cell(current_batch_size, 1); % [cell] Temporarily holds the assigned 3D voxel data for each block.
    batch_indices = cell(current_batch_size, 1); % [cell] Stores the physical matrix indices where the block should be re-inserted later.
	
    batch_grain_counts = zeros(current_batch_size, 1); % [array] Tracks how many grains were evaluated per block for metrics.

    % Extract grid parameters so parfor loop has clean, read-only variables to prevent scope issues
    xvec_local = xvec;
    yvec_local = yvec;
    zvec_local = zvec;
    x_edges_local = x_edges;
    y_edges_local = y_edges;
    z_edges_local = z_edges;

    % Parallel loop over blocks in batch (Actual Grain Assignment)
    parfor local_idx = 1:current_batch_size
        bx = block_list_batch(local_idx, 1);
        by = block_list_batch(local_idx, 2);
        bz = block_list_batch(local_idx, 3);

        % Compute index ranges for voxels within the current block
        ix1 = (bx-1)*block_size + 1; ix2 = min(bx*block_size, Nx);
        iy1 = (by-1)*block_size + 1; iy2 = min(by*block_size, Ny);
        iz1 = (bz-1)*block_size + 1; iz2 = min(bz*block_size, Nz);

        % Block world coordinates (physical coordinates at block boundaries)
        block_xmin = xvec_local(ix1); block_xmax = xvec_local(ix2);
        block_ymin = yvec_local(iy1); block_ymax = yvec_local(iy2);
        block_zmin = zvec_local(iz1); block_zmax = zvec_local(iz2);

        % Compute relevant grains directly without accessing spatial_grid
        relevant_grains = find(...
            grain_seeds_single(:,1) >= (block_xmin - max_influence) & ...
            grain_seeds_single(:,1) <= (block_xmax + max_influence) & ...
            grain_seeds_single(:,2) >= (block_ymin - max_influence) & ...
            grain_seeds_single(:,2) <= (block_ymax + max_influence) & ...
            grain_seeds_single(:,3) >= (block_zmin - max_influence) & ...
            grain_seeds_single(:,3) <= (block_zmax + max_influence)); % [array] List of grains near enough to affect voxels in THIS specific block.


        % If block contains any relevant grains, proceed with assignment
        if ~isempty(relevant_grains)
		
            % Extract actual physical 3D coordinates for every voxel in this block
            [Xb, Yb, Zb] = ndgrid(single(xvec_local(ix1:ix2)), single(yvec_local(iy1:iy2)), single(zvec_local(iz1:iz2)));

            % Initialize arrays to hold the assigned Grain IDs and calculated distances
            block_phase = zeros(size(Xb), 'uint32'); % [3D matrix] Local map of Grain IDs for this block.
			
            min_dist2 = inf(size(Xb), 'single'); % [3D matrix] Tracks the closest grain found so far for every voxel (starts at infinity).

            % For each relevant grain, assign voxels in block to their nearest grain
            for k = relevant_grains'
			
                % Calculate delta distance from voxel coordinate to grain center			
                dX = Xb - grain_seeds_single(k,1);
                dY = Yb - grain_seeds_single(k,2);
                dZ = Zb - grain_seeds_single(k,3);

                % Compute an ellipsoidal distance metric (factors in aspect ratio)
                current_dist2 = (dX.^2) * inv_r_xy2(k) + (dY.^2) * inv_r_xy2(k) + (dZ.^2) * inv_r_z2(k); % [3D matrix] Distance from current grain to every voxel in the block.

                % If the current grain is closer than any previously checked grain, update the mask
                update_mask = current_dist2 < min_dist2; % [3D logical] True where the current grain 'k' is the closest found so far.

                % Assign the current Grain ID (k) to voxels where it won
                block_phase(update_mask) = k;
				
                % Update the minimum distance tracker				
                min_dist2(update_mask) = current_dist2(update_mask);
            end

            % Save completed block data back to cell arrays
            batch_results{local_idx} = block_phase;
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = length(relevant_grains);
        else
            % If no grains were near this block, leave it empty		
            batch_results{local_idx} = [];
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = 0;
        end
    end

    % --- Data Reassembly ---
    % All blocks in the current batch have been processed in parallel. 
    % Insert the small blocks back into the main phase cube sequentially.
    batch_grains_processed = 0;
    batch_blocks_with_grains = 0;

    for local_idx = 1:current_batch_size
        if ~isempty(batch_results{local_idx})

            % Retrieve coordinates where this block belongs
            indices = batch_indices{local_idx};
            ix1 = indices(1); ix2 = indices(2);
            iy1 = indices(3); iy2 = indices(4);
            iz1 = indices(5); iz2 = indices(6);

            % Splice the block data into the main 3D matrix
            phase(ix1:ix2, iy1:iy2, iz1:iz2) = cast(batch_results{local_idx}, phase_data_type);
			
            batch_blocks_with_grains = batch_blocks_with_grains + 1;
        end
        % Sum up grains processed in batch for metrics
        batch_grains_processed = batch_grains_processed + batch_grain_counts(local_idx);
    end

    % Update global progress summary counters 
    completed_blocks = completed_blocks + current_batch_size;
    total_grains_processed = total_grains_processed + batch_grains_processed;
    blocks_with_grains = blocks_with_grains + batch_blocks_with_grains;

    % Record time spent on this batch to estimate future progress
    batch_time = toc(batch_timer);
    batch_times(end+1) = batch_time;

    % Calculate overall statistics for progress reporting
    progress_pct = 100 * completed_blocks / total_blocks;
    elapsed_total = toc(main_timer);
    avg_grains_per_block = batch_grains_processed / current_batch_size;

    % Estimate Time of Arrival (ETA) for remaining batches based on average speed
    if batch_idx > 1
        avg_batch_time = mean(batch_times);
        remaining_batches = n_batches - batch_idx;
        eta_seconds = remaining_batches * avg_batch_time;
    else
        eta_seconds = batch_time * (n_batches - 1);
    end

    % Print simple progress line to console
    fprintf('Batch %d/%d: Blocks %d-%d (%.1f%%) | %.1f grains/block avg | Elapsed: %.1fs | ETA: %.1fs (%.1f min)\n', ...
        batch_idx, n_batches, batch_start, batch_end, progress_pct, avg_grains_per_block, elapsed_total, eta_seconds, eta_seconds/60);

    % Every few batches, print a detailed summary of system performance
    if mod(batch_idx, max(1, floor(n_batches/5))) == 0 || batch_idx == n_batches
	
        overall_avg_grains = total_grains_processed / completed_blocks;
        processing_rate = completed_blocks / elapsed_total;

        % Speedup factor calculates efficiency gained by not checking every single grain for every block		
        speedup_factor = N_grains / overall_avg_grains;

        fprintf('  --> Detailed Stats: %.1f grains/block overall | %.1f blocks/sec | %.1fx speedup from spatial filtering\n', ...
            overall_avg_grains, processing_rate, speedup_factor);

        if batch_idx < n_batches
            fprintf('  --> Memory: %d/%d blocks processed | %d/%d non-empty | %.1f%% efficiency\n', ...
                completed_blocks, total_blocks, blocks_with_grains, completed_blocks, 100*blocks_with_grains/completed_blocks);
        end
    end
end

% 5-5 Final Performance Summary
total_time = toc(main_timer);
avg_grains_per_block = total_grains_processed / completed_blocks;
speedup_factor = N_grains / avg_grains_per_block;
processing_rate = double(Nx*Ny*Nz) / total_time;

fprintf('\n=== COMPLETION SUMMARY ===\n');
fprintf('Total time: %.1f seconds (%.2f minutes)\n', total_time, total_time/60);
fprintf('Processing rate: %s voxels/second\n', regexprep(sprintf('%.0f', processing_rate),'(\d)(?=(\d{3})+$)','$1,'));
fprintf('\nSpatial filtering efficiency:\n');
fprintf('  Average grains per block: %.1f (vs %d total grains)\n', avg_grains_per_block, N_grains);
fprintf('  Computational speedup: %.1fx\n', speedup_factor);
fprintf('  Blocks processed: %d/%d\n', completed_blocks, total_blocks);
fprintf('  Non-empty blocks: %d (%.1f%%)\n', blocks_with_grains, 100*blocks_with_grains/completed_blocks);

fprintf('\nParallel processing efficiency:\n');
fprintf('  Workers used: %d\n', num_workers);
fprintf('  Batches processed: %d\n', n_batches);
fprintf('  Average batch time: %.1f seconds\n', mean(batch_times));

% Verify that all voxels have grain assignments; warn if any left blank
unassigned = sum(phase(:) == 0);
if unassigned > 0
    fprintf('\nWARNING: %d unassigned voxels detected!\n', unassigned);
else
    fprintf('\nSUCCESS: All %s voxels assigned successfully!\n', ...
        regexprep(sprintf('%.0f',double(Nx*Ny*Nz)),'(\d)(?=(\d{3})+$)','$1,'));
end

fprintf('\nParallel + Spatial filtering microstructure generation completed!\n');



%% 6. Overwrite With Pore Regions -----------------------------------------

% Setup & Pre-calculation
sz_phase = size(phase); % [array] 1x3 vector containing the dimensions of the phase matrix [Nx, Ny, Nz].
Nx = sz_phase(1); Ny = sz_phase(2); Nz = sz_phase(3);

% Prepare storage for parallel execution
pore_linear_indices = cell(N_pores, 1); % [cell array] Stores the 1D global array indices of all voxels that belong to each pore.

% Only proceed if pore positions exist
if ~isempty(x_pore)
    fprintf('Inserting %d Keyhole Voids (Vertical Prolate Spheroids)...\n', N_pores);

    for k = 1:N_pores
        % --- A. Geometry: Penny Shape ---
        c = D_Feret_max(k) / 2; % Half-axis in Z (Height)

        if c <= 0 || V_pore(k) <= 0
            continue;
        end

        % Calculate Radius (a = b) required to match the experimental Volume.
        % Formula: Volume = (4/3) * pi * a * b * c
        % Since a = b:  a = sqrt( Volume / ((4/3) * pi * c) )
        val_sq = V_pore(k) / ((4/3) * pi * c);

        if val_sq <= 0
            continue;
        end

        a = sqrt(val_sq);
        b = a;

        % --- B. Physics Enforcement ---
        % Keyhole pores are deep and narrow (c > a).
        % If the derived width (a) is larger than the height (c),
        % the data likely represents a gas pore that was misclassified.
        if a > c
            r_eq = ((3 * V_pore(k)) / (4 * pi))^(1/3);
            a = r_eq;
            b = r_eq;
            c = r_eq;
        end

        % --- C. Bounding Box ---
        % Calculate the maximum possible radius needed to enclose the rotated shape.
        buffer = 1.1;
        max_r = max([a, b, c]) * buffer;

        % Get indices within the bounding box
        idx_x = find(xvec >= (x_pore(k) - max_r) & xvec <= (x_pore(k) + max_r));
        idx_y = find(yvec >= (y_pore(k) - max_r) & yvec <= (y_pore(k) + max_r));
        idx_z = find(zvec >= (z_pore(k) - max_r) & zvec <= (z_pore(k) + max_r));

        if isempty(idx_x) || isempty(idx_y) || isempty(idx_z)
            continue;
        end
		
        [X_sub, Y_sub, Z_sub] = ndgrid(xvec(idx_x), yvec(idx_y), zvec(idx_z));
        Xc = X_sub - x_pore(k);
        Yc = Y_sub - y_pore(k);
        Zc = Z_sub - z_pore(k);

        % --- D. Ellipsoid Evaluation ---
        % Implicit equation: (x/a)^2 + (y/b)^2 + (z/c)^2 <= 1
        mask_sub = (Xc/a).^2 + (Yc/b).^2 + (Zc/c).^2 <= 1;

        % --- E. Index Mapping ---
        % Find subscripts relative to the small sub-grid
        [ii, jj, kk] = ind2sub(size(mask_sub), find(mask_sub));

        % Map to Global Indices
        if ~isempty(ii)
            % Map sub-grid subscripts to Global Grid indices
            global_x = idx_x(ii);
            global_y = idx_y(jj);
            global_z = idx_z(kk);

            % Store Linear Indices
            pore_linear_indices{k} = sub2ind(sz_phase, ...
                double(global_x(:)), ...
                double(global_y(:)), ...
                double(global_z(:)));
        end

    end

    all_inds = vertcat(pore_linear_indices{:});
    % Assign the Keyhole Phase ID
    phase(all_inds) = N_grains + 1;

    % Statistics
    active_pores = sum(~cellfun(@isempty, pore_linear_indices));
    fprintf('Keyhole Insertion Complete. %d active pores inserted.\n', active_pores);
    nphase = N_grains+1;
else
    nphase = N_grains;
    fprintf('No Keyhole data found. Skipping insertion.\n');
end

fprintf('Phase assignment completed.\n');

%% 7. Assign Orientations -------------------------------------------------

ori_vec = zeros(3*N_grains, 1); % [array] A 1D array of size (3*N_grains) x 1 to store the 3-component Rodrigues vector for each grain sequentially.


% If assigning a fiber texture (preferred orientation about a direction)
if strcmp(orientation_type, 'textured')
    % --- CASE  EOS/Keyhole (Fiber Texture) ---
    % Simulates columnar growth: Crystal [101] || Sample Z
    fprintf('Generating [101] Fiber Texture \n');

    % Target: Align Crystal [101] to Sample Z [001]
    c_dir = [1 0 1]; c_dir = c_dir / norm(c_dir); % [array] The chosen principal axis of the crystal lattice.
    s_dir = [0 0 1];                              % [array] The physical build direction of the sample.

    % Base Rotation (Aligns Crystal [101] to Sample Z)
    % Calculate axis-angle representation for the shortest rotation between c_dir and s_dir using cross and dot products.	
    v = cross(c_dir, s_dir);
    s = norm(v);
    c = dot(c_dir, s_dir);
	
    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]; % [matrix] Skew-symmetric cross-product matrix of v.
	
    if s < 1e-8
        % If vectors are already perfectly aligned, use the Identity matrix.	
        R_base = eye(3);
    else
        % Rodrigues' rotation formula to generate the base rotation matrix.	    
        R_base = eye(3) + vx + vx^2*((1-c)/s^2);  % [3x3 matrix] Matrix that explicitly rotates [101] into [001].
    end

    for k = 1:N_grains
        % Random Spin (0-360) around Sample Z
        % This creates the "Ring" in the Pole Figure (Fiber Texture)
        phi = rand() * 2*pi;
		
        R_spin = [cos(phi) -sin(phi) 0;
            sin(phi)  cos(phi) 0;
            0         0       1]; % [3x3 matrix] Rotation matrix for the random azimuthal spin.

        % Gaussian Wobble (Spread)
        % Random axis in XY plane, small angle
        wobble_ang = deg2rad(sigma_spread * randn());
        wb_ax_ang  = rand() * 2*pi;
        ax_w = [cos(wb_ax_ang); sin(wb_ax_ang); 0];
		
        K = [0 -ax_w(3) ax_w(2); ax_w(3) 0 -ax_w(1); -ax_w(2) ax_w(1) 0];
		
        R_wobble = eye(3) + sin(wobble_ang)*K + (1-cos(wobble_ang))*K*K; % [3x3 matrix] Small perturbation matrix to simulate non-perfect epitaxial alignment.

        % Combine: Spin -> Wobble -> Base Alignment
        % Note: Matrices multiply Right-to-Left in standard Column vector notation
        % x_new = R_spin * R_wobble * R_base * x_old
        R_total = R_spin * R_wobble * R_base; % [3x3 matrix] The final combined rotation matrix for this specific grain.

        % Convert to Rodrigues
        tr = trace(R_total);
        theta = acos(max(min((tr-1)/2, 1), -1));
		
        if abs(theta) < 1e-6
            rod = [0 0 0];
        else
            r_axis = (1/(2*sin(theta))) * ...
                [R_total(3,2)-R_total(2,3);
                R_total(1,3)-R_total(3,1);
                R_total(2,1)-R_total(1,2)];
            rod = r_axis * tan(theta/2);
        end
		
        % Store the resulting 3-component vector in the flat array		
        idx = 3*(k-1)+1;
        ori_vec(idx:idx+2) = rod(:);
    end

elseif strcmp(orientation_type, 'random')
    % --- CASE : Lack of Fusion (Random Texture) ---
    % Simulates unmelted powder or disrupted epitaxy
    % Uses Uniform Random Rotation (Haar Measure)

    fprintf('Generating Random Texture \n');
    for k = 1:N_grains
        % Generate Uniform Random Quaternion
		% Generate three independent uniform random variables between 0 and 1
        u1 = rand(); u2 = rand(); u3 = rand();

        % Calculate quaternion components mapping uniform variables to spherical space
        q_w = sqrt(1-u1) * sin(2*pi*u2);
        q_x = sqrt(1-u1) * cos(2*pi*u2);
        q_y = sqrt(u1)   * sin(2*pi*u3);
        q_z = sqrt(u1)   * cos(2*pi*u3);

        % Convert to Rodrigues
        % Rodrigues vector = [qx, qy, qz] / qw		
        if abs(q_w) < 1e-8
            rod_vec = [q_x q_y q_z] / 1e-9; % Avoid Division by Zero (Inf)
        else
            rod_vec = [q_x q_y q_z] / q_w;
        end

        % Store the resulting 3-component vector in the flat array
        idx = 3*(k-1)+1;
        ori_vec(idx:idx+2) = rod_vec;
    end
end




%% 8. --- Preparation for IPF coloring ---

if strcmp(orientation_type, 'textured')

    cs = crystalSymmetry('m-3m'); % [MTEX object] Defines the crystal symmetry as cubic (Laue group m-3m).
	
    rods = reshape(ori_vec(1:3*N_grains), [3, N_grains])'; % [matrix] Reshapes the flat 1D ori_vec into an N_grains x 3 matrix.
	
    ori_mtex = orientation.byRodrigues(rods, cs); % [MTEX object] Converts the numeric Rodrigues vectors into MTEX orientation objects.

     % --- Generate IPF Colors ---
    RD = vector3d.Z;  % [MTEX object] Reference direction for the IPF color mapping (Sample Z-axis).
	
    ipf_key = ipfHSVKey(cs); % [MTEX object] Initializes the HSV color key for cubic symmetry.
    ipf_key.inversePoleFigureDirection = RD; % Set Z as the reference
	
    cmap_ipf = ipf_key.orientation2color(ori_mtex); % [matrix] N_grains x 3 array containing the RGB color triplet for every grain.


elseif strcmp(orientation_type, 'random')

    cs = crystalSymmetry('m-3m');
    rods = reshape(ori_vec(1:3*N_grains), [3, N_grains])';
    ori_mtex = orientation('rodrigues', rods, cs);
    ipf_key = ipfColorKey(cs);
    ipf_key.inversePoleFigureDirection = vector3d.Z;
    cmap_ipf = ipf_key.orientation2color(ori_mtex); % Generate IPF colors
end


%% 9. Save To HDF5 --------------------------------------------------------

n_phases = max(phase(:)); % [integer] The absolute maximum phase ID found in the volume (Total Grains + 1 Void Phase).

% --- Save to MATLAB format ---
if strcmp(save_data, 'on')

    % Save only the essential grid parameters and basic geometry stats for quick reloading
    save('xct_poly_params_SS316L.mat','xvec','yvec','zvec','Nx','Ny','Nz', ...
        'cube_len','N_grains','n_phases');
		
    % Save all current workspace variables to a .mat file (version 7.3 for large arrays)
    save('xct_poly_params_SS316L_full.mat','-v7.3')
	
elseif strcmp(save_data, 'off')
    fprintf('save_data off');
end

% --- Save to HDF5 format ---
if strcmp(elas3d_input, 'on')

    h5filename = 'input_structure_poly.h5';
    % Delete previous file if it exists (overwrite protection)
    if isfile(h5filename)
        delete(h5filename);
        fprintf('Previous file "%s" deleted.\n', h5filename);
    end

    % Flatten phase array into a vector for writing into HDF5
    pix_flat = int32(phase(:));
    % Create HDF5 dataset '/pix' to store the 1D phase information
    h5create(h5filename, '/pix', [length(pix_flat) 1], 'Datatype', 'int32');
    % Write the flattened phase array into the HDF5 file
    h5write(h5filename, '/pix', pix_flat);
    % Annotate '/pix' dataset with an attribute containing its original 3D dimensions [Nx, Ny, Nz]
    h5writeatt(h5filename, '/pix', 'dimensions', int32([Nx Ny Nz]));
    % Create HDF5 dataset '/orientation' for the 1D list of Rodrigues vectors
    h5create(h5filename, '/orientation', size(ori_vec), 'Datatype', 'double');
    % Write orientation data to the file
    h5write(h5filename, '/orientation', ori_vec);
    % Annotate the '/orientation' dataset with its array dimensions
    h5writeatt(h5filename, '/orientation', 'dimensions', int32(size(ori_vec)));

    disp(['Saved HDF5 input structure to ', h5filename]);

elseif strcmp(elas3d_input, 'off')
    fprintf('HDF5 OFF \n');
end


%% 10. Output Summary  --------------------------------------------------------

% Open file for writing
fid = fopen('output_summary.txt', 'w');

% Print summary of grains and pores found
fprintf(fid, 'Total number of Grains: %s, Pores: %s\n', ...
    regexprep(sprintf('%.0f',N_grains),'(\d)(?=(\d{3})+$)','$1,'), ...
    regexprep(sprintf('%.0f',nphase-N_grains),'(\d)(?=(\d{3})+$)','$1,'));

% Print physical size of the simulation domain and voxel resolution
fprintf(fid, 'Domain size: %.1f x %.1f x %.1f um^3\n', ...
    dx*1e3, dx*1e3, dx*1e3);
fprintf(fid, 'Domain grid: %d x %d x %d = %s voxels\n', Nx, Ny, Nz, ...
    regexprep(sprintf('%.0f', Nx*Ny*Nz), '(\d)(?=(\d{3})+$)', '$1,'));
fprintf(fid, 'Voxel size: %.5f mm (%.2f µm)\n', dx, dx*1e3);

% --- Voxel Counting (Volume Fractions / Porosity) ---
% Count the number of voxels belonging to grains
num_vox_grain = sum( phase(:) >= 1 & phase(:) <= N_grains ); % [integer] Total voxels assigned a valid grain ID (1 to N_grains).

% Count the number of voxels assigned as pores
num_vox_pore = sum( phase(:) == (N_grains+1) ); % [integer] Total voxels assigned the unique void ID (N_grains + 1).

% Total voxel count
num_vox_total = numel(phase); % [integer] Absolute number of elements in the phase matrix.

% Print percentages for grains vs. pores
fprintf(fid, 'Grain voxels: %d (%.3f%%)\n', num_vox_grain, 100*num_vox_grain/num_vox_total);
fprintf(fid, 'Pore  voxels: %d (%.3f%%)\n', num_vox_pore, 100*num_vox_pore/num_vox_total);
fprintf(fid, 'Total voxels: %d\n', num_vox_total);

% Print mean grain diameter (in microns, for clarity)
fprintf(fid, 'Mean grain size (XY): %.2f um\n', mean(grain_diameters)*1e3);

% Print mean aspect ratio (Z/XY)
fprintf(fid, 'Mean aspect ratio (Z/XY): %.2f\n', mean(aspect_ratios));

% Final completion message
fprintf(fid, 'Done! Phase voxel mesh saved for full cuboid domain with defects.\n');

% Close file
fclose(fid);

disp('Done!');
tElapsed = toc(tStart);

% --- Performance Profiling ---
% Calculate the total elapsed time since 'tic' was called at 'tStart'
fprintf('\nComputation completed in:\n');
fprintf('  %.2f seconds\n', tElapsed);
fprintf('  %.2f minutes\n', tElapsed/60);
fprintf('  %.2f hours\n', tElapsed/3600);

%% 11. Plot --------------------------------------------------------------

if strcmp(micro_plot, 'on')

    if strcmp(microstructure3d_type, 'surface')
        % --- Plot Outer Shell Only (Faster Rendering) ---
        [Nx,Ny,Nz] = size(phase);
		
        % Create a logical mask identifying only the 6 exterior faces of the cube
        surface_mask = false(Nx,Ny,Nz);
        surface_mask(1,:,:)   = true; surface_mask(Nx,:,:)  = true;
        surface_mask(:,1,:)   = true; surface_mask(:,Ny,:)  = true;
        surface_mask(:,:,1)   = true; surface_mask(:,:,Nz)  = true;
		
        % Flatten 3D coordinate arrays to 1D for the scatter3 function
        phase_flat = phase(:);
        X_flat = X(:); Y_flat = Y(:); Z_flat = Z(:);
		
        % Extract indices for only the surface voxels
        surface_inds = find(surface_mask(:));
		
        % Apply stride subsampling to further speed up plotting
        surface_inds = surface_inds(1:stride:end);
        grain_inds = phase_flat(surface_inds);

        % Filter out unassigned space and pore voxels
        valid = grain_inds >= 1 & grain_inds <= size(cmap_ipf,1);
        surface_inds = surface_inds(valid);
        grain_inds = grain_inds(valid);

        % Map each point to its corresponding IPF RGB color
        pt_colors = cmap_ipf(grain_inds, :);

        figure; clf; hold on;
        fprintf('Plotting 3D microstructure SURFACE by scatter3 (IPF color)...\n');

        scatter3(X_flat(surface_inds), Y_flat(surface_inds), Z_flat(surface_inds), 9, ...
            pt_colors, 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0);

        % Overlay Pores as distinct spheres
        if exist('n_pores','var') && n_pores > 0
            [xx,yy,zz]=sphere(30);
            for i=1:n_pores
                surf(x(i)+r(i)*xx, y(i)+r(i)*yy, z(i)+r(i)*zz, ...
                    'FaceAlpha',0.7,'EdgeColor','none','FaceColor',[0.4,0.4,0.6]);
            end
        end

        hx = xlabel('X_1 (mm)'); set(hx, 'Rotation', 21);
        hy = ylabel('X_2 (mm)'); set(hy, 'Rotation', -21);
        zlabel('X_3 (mm)');
        title('3D Grain Microstructure (IPF color FCC [101] RD)');
        view([-42.71 22.97]);
        axis equal tight; grid on; box on;

        xticks(linspace(xmin, xmax, 6));
        yticks(linspace(ymin, ymax, 6));
        zticks(linspace(zmin, zmax, 6));
        xlim([xmin,xmax]);
        ylim([ymin,ymax]);
        zlim([zmin,zmax]);
        set(gca,'FontSize',15,'LineWidth',1.2);

        hAxes = findobj(gcf,"Type","axes");
        hAxes.BoxStyle = "full";
        hAxes.ClippingStyle = "3dbox";
        hold off

    elseif strcmp(microstructure3d_type, 'full3d')
        % --- Plot Full Dense 3D Volume ---
        figure; clf; hold on;
        fprintf('Plotting 3D microstructure by scatter3 (IPF color)...\n');
		
        for k = 1:N_grains
            inds = find(phase(:)==k);
            if isempty(inds), continue; end
            inds = inds(1:stride:end);
			
            % Plot all internal points for this specific grain			
            scatter3(X(inds), Y(inds), Z(inds), 9, ...
                repmat(cmap_ipf(k,:), numel(inds),1), ...
                'filled', 'MarkerFaceAlpha',0.6, 'MarkerEdgeAlpha',0.05);
        end

        % Plot pores
        if exist('n_pores','var') && n_pores > 0
            [xx,yy,zz]=sphere(30);
            for i=1:n_pores
                surf(x(i)+r(i)*xx, y(i)+r(i)*yy, z(i)+r(i)*zz, ...
                    'FaceAlpha',0.7,'EdgeColor','none','FaceColor',[0.4,0.4,0.6]);
            end
        end

        hx = xlabel('X_1 (mm)'); set(hx, 'Rotation', 21);
        hy = ylabel('X_2 (mm)'); set(hy, 'Rotation', -21);
        zlabel('X_3 (mm)');
        xlim([xmin cube_len]); ylim([ymin cube_len]); zlim([zmin cube_len]);
        title('3D Grain Microstructure (IPF color FCC [101] RD)');
        view([-42.71 22.97]);
        axis equal tight; grid on; box on;

        xticks(linspace(xmin, xmax, 6));
        yticks(linspace(ymin, ymax, 6));
        zticks(linspace(zmin, zmax, 6));
        xlim([xmin,xmax]);
        ylim([ymin,ymax]);
        zlim([zmin,zmax]);

        set(gca,'FontSize',15,'LineWidth',1.2);
        hAxes = findobj(gcf,"Type","axes");
        hAxes.BoxStyle = "full";
        hAxes.ClippingStyle = "3dbox";
        hold off
		
    end
	

    %% --- Crystallographic Pole Figures (MTEX) ---

    % Randomly sample a subset (K) of grain orientations to compute the pole figure efficiently
    N = numel(ori_vec) / 3; 
    K = min(K_grain, N);
    idx_grains = randperm(N, K);
	
    % Extract Rodrigues vectors for the sampled grains
    ori_idx = bsxfun(@plus, 3*(idx_grains'-1), (1:3));
    ori_idx = ori_idx';
    ori_idx = ori_idx(:);
    ori_vec_downsampled = ori_vec(ori_idx);
    ori_vec_matrix = reshape(ori_vec_downsampled, 3, []).';
	
    % Convert to MTEX orientation objects
    ori = orientation.byRodrigues(ori_vec_matrix, cs, ss);
	
    % Estimate Orientation Distribution Function (ODF) using a von Mises-Fisher kernel
    psi = SO3vonMisesFisherKernel('halfwidth',5*degree);
    % Estimate the unimodal ODF from the sampled orientations using the kernel
    odf = unimodalODF(ori, psi);
	
    % Set coordinate conventions for pole figure display
    setMTEXpref('xAxisDirection','east');
    setMTEXpref('zAxisDirection','outOfPlane');

    % --- Plot Standard Pole Figure (PF) ---
    figure; clf;
    plotPDF(odf, [Miller(1,0,0,cs), Miller(1,1,0,cs), Miller(1,1,1,cs)], ...
        'antipodal', 'contourf', 'smooth', 'resolution', 1*degree, 'minmax');
    colorbar;


elseif strcmp(micro_plot, 'off')
    disp('micro_plot_disabled')
end

%% ---- END ----