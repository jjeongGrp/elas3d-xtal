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
format compact; clear; close all;
tStart = tic; % Start a timer to track total microstructure generation time

%% 1. Configuration & Parameters -------------------------------------------------------------
% ==============================================================================
% Defines the physical and computational dimensions of the domain, the
% statistical properties of the grains (size and aspect ratio), crystallographic
% texture settings, and flags controlling outputs and visualizations.
% ==============================================================================

% --- Domain Variables ---
max_len = 5.0;           % [mm] max_len: Absolute maximum physical length allowed for the domain boundaries.
cube_len = 5.0;          % [mm] cube_len: The actual physical side length of the cubic domain being generated.
grid_num = 500;          % [voxels] grid_num: Number of voxels along each 1D edge of the cube (Resolution).
dx = cube_len/grid_num;  % [mm] dx: The physical length of a single voxel edge.
xmin = 0;                % [mm] xmin: Lower boundary of the domain on the X-axis (Origin).
ymin = 0;                % [mm] ymin: Lower boundary of the domain on the Y-axis (Origin).
zmin = 0;                % [mm] zmin: Lower boundary of the domain on the Z-axis (Origin).

% --- Grain Statistics Variables ---
mean_d = 0.055;        % [mm] mean_d: Target average equivalent diameter of the grains in the XY plane.
std_d  = 0.35*mean_d;  % [mm] std_d: Standard deviation of the grain diameter (controls size variability).


% Aspect ratio variables (defines if grains are equiaxed (~1.0) or columnar (>1.0))
mean_aspect = 3.18;              % [ratio] mean_aspect: Average length-to-width ratio of the grains.
std_aspect = 0.25 * mean_aspect; % [ratio] std_aspect: Standard deviation of the aspect ratio.


% --- Orientation Settings Variables ---
orientation_type = 'textured'; % [string] orientation_type: 'textured' for aligned orientations, or 'random'.
sigma_spread = 15;             % [degrees] sigma_spread: If 'textured', defines how much orientations can deviate from the primary axis.

% --- General Descriptor Variables ---
pore_type = 'EOS';             % [string] pore_type: Label to identify the source/type of pore data.
desc = 'Spherical Gas Pore';   % [string] desc: Description of the defect morphology.

% --- Output Settings Variables ---
save_data = 'on';              % [string] save_data: Flag to save MATLAB workspace variables to a .mat file. 'on' or 'off'
elas3d_input = 'on';           % [string] elas3d_input: Flag to generate and save the .h5 file for ELAS3D. 'on' or 'off'
phase_data_type = 'uint32';    % [string] phase_data_type: Data type for the 3D matrix storing grain IDs (uint32 supports up to ~4.2 billion grains)


% --- Visualization Variables ---
stride = 1;                        % [integer] stride: Downsampling factor for scatter plots (higher = faster plotting but lower visual resolution).
microstructure3d_type = 'surface'; % [string] microstructure3d_type: 'surface' plots outer faces only, 'full3d' plots internal volumes.
max_pores_show = 1;    % [integer] max_pores_show: Limits how many individual pore isosurfaces are rendered.
num_grains_iso = 10;   % [integer] num_grains_iso: Limits how many individual grain isosurfaces are rendered.
micro_plot = 'on';     % [string] micro_plot: Master flag to turn 3D rendering 'on' or 'off'.
K_grain = 2000;        % [integer] number of grains to randomly sample for pole figure analysis, limited to available grains.

% --- Pore Settings Variables ---
pore = 'off';          % [string] pore: Master flag to include or exclude defects in the generated volume. 'on' or 'off'
input_pore_file = 'Defect_Locations_EOS.xlsx'; % [string] input_pore_file: Filename containing experimental/synthetic defect coordinates.

%% 2. Pore Data -----------------------------------------------------------
% ==============================================================================
% Reads explicit pore coordinates and size geometries from an external spreadsheet.
% If the file is missing, it skips pore generation and initializes empty arrays
% to prevent the script from crashing downstream.
% ==============================================================================

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
% ==============================================================================
% Constructs the 3D spatial grid (voxel mesh) for the simulation domain.
% It defines the physical boundaries, generates the 1D coordinate vectors,
% and uses ndgrid to create the full 3D coordinate matrices. 
% It also includes (currently commented out) diagnostic checks to verify grid 
% extents and calculates how many experimental pores are large enough to be 
% resolved by the chosen voxel size.
% ==============================================================================

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
% ==============================================================================
% Calculates the required number of grains to fill the domain based on 
% target volume estimates. It then uses lognormal statistical distributions 
% to assign a unique size (XY diameter) and aspect ratio (Z-axis elongation) 
% to every single grain. Finally, it randomly scatters the seed coordinates
% (centers) of these grains throughout the 3D space.
% ==============================================================================

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
% ==============================================================================
% Begins the core voxel-to-grain assignment process. Because comparing 
% millions of voxels to thousands of grains is computationally expensive, 
% this section implements two major optimizations: 
% 1. Multi-threading: It spins up a parallel pool using all available CPU cores.
% 2. Spatial Acceleration: It creates a coarse 'macro-grid'. Instead of a 
% voxel checking every grain in the domain, it only checks grains that 
% physically overlap with its local macro-grid cell.
% It then breaks the domain into smaller 3D blocks and distributes them across
% parallel workers (parfor loops) to assign each voxel to its closest grain.
% ==============================================================================

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
% ==============================================================================
% Divides the entire domain into a coarse 15x15x15 grid. For each coarse cell, 
% it compiles a short list of "relevant grains" whose boundaries are close 
% enough to potentially touch voxels inside that cell.
% ==============================================================================
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
% ==============================================================================
% Breaks the main, massive voxel grid into smaller, manageable 3D chunks. 
% These blocks are then grouped into "batches" so they can be fed to 
% the parallel workers efficiently without overwhelming system memory.
% ==============================================================================

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
% ==============================================================================
% The main parallel execution loop. It iterates through the batches, assigns 
% blocks to parallel workers, determines which grains are relevant to that 
% specific block, calculates distances, and assigns the final Grain ID to 
% the voxels. Then, it merges the completed block back into the main matrix.
% ==============================================================================
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
% ==============================================================================
% Takes the fully dense grain structure generated in the previous steps 
% and introduces defects. To do this efficiently, it iterates through 
% the experimental pore coordinates and uses a bounding box optimization:
% instead of checking every voxel in the domain to see if it 
% falls inside a pore, it only checks the small neighborhood of voxels 
% immediately surrounding each pore's center. Voxels falling inside the pore 
% radius are overwritten with a unique ID representing void space.
% ==============================================================================


% Setup & Pre-calculation
sz_phase = size(phase); % [array] 1x3 vector containing the dimensions of the phase matrix [Nx, Ny, Nz].
Nx = sz_phase(1); Ny = sz_phase(2); Nz = sz_phase(3);

% Prepare storage for parallel execution
pore_linear_indices = cell(N_pores, 1); % [cell array] Stores the 1D global array indices of all voxels that belong to each pore.

% Only proceed if pore positions exist
if ~isempty(x_pore)
    fprintf('Inserting %d Gas Voids (Spherical)... \n', N_pores);

    % Iterate through every pore from the input dataset
    for k = 1:N_pores

        % --- Geometry: Sphere ---
        r_pore = d(k) / 2.0; % [float] Radius of the current pore [mm].

        % Skip invalid or zero-size pores
        if r_pore <= 0 || isnan(r_pore)
            continue;
        end

        % --- Bounding Box ---
        % Only search voxels within +/- r of the center
        buffer = 1.1;         % [float] Multiplier to slightly expand the search box to ensure voxel edges are caught.
		
        max_r = r_pore * buffer; % [float] The expanded search radius [mm].

        % Calculate physical boundaries of the local search box
        x_min = x_pore(k) - max_r;
        x_max = x_pore(k) + max_r;
        y_min = y_pore(k) - max_r;
        y_max = y_pore(k) + max_r;
        z_min = z_pore(k) - max_r;
        z_max = z_pore(k) + max_r;

        % Find the vector indices of the grid points that fall inside the bounding box
        idx_x = find(xvec >= x_min & xvec <= x_max);
        idx_y = find(yvec >= y_min & yvec <= y_max);
        idx_z = find(zvec >= z_min & zvec <= z_max);

        % If the bounding box is completely outside the simulated domain, skip this pore
        if isempty(idx_x) || isempty(idx_y) || isempty(idx_z)
            continue;
        end

        % Create local grid
        [X_sub, Y_sub, Z_sub] = ndgrid(xvec(idx_x), yvec(idx_y), zvec(idx_z));         % [3D matrices] Local 3D coordinate grids containing only the voxels near the pore.

        % Shift the local coordinates so the pore center is at (0,0,0)
        Xc = X_sub - x_pore(k);
        Yc = Y_sub - y_pore(k);
        Zc = Z_sub - z_pore(k);

        % --- Sphere Equation ---
        mask_sub = (Xc.^2 + Yc.^2 + Zc.^2) <= r_pore^2; % [3D logical] True (1) if the voxel center is physically inside the spherical pore volume.

        % Convert the logical mask into local 3D matrix subscripts
        [ii, jj, kk] = ind2sub(size(mask_sub), find(mask_sub));

        % Map to Global Indices
        if ~isempty(ii)
		
            % Convert local subscripts back into the global 1D vector indices
            global_x = idx_x(ii);
            global_y = idx_y(jj);
            global_z = idx_z(kk);

            % --- Strict Column Vectors ---
            % Force all inputs to be column vectors using (:).
            pore_linear_indices{k} = sub2ind(sz_phase, ...
                double(global_x(:)), ...
                double(global_y(:)), ...
                double(global_z(:))); % [array] Converts global 3D subscripts (X,Y,Z) into a single 1D linear index for fast matrix assignment.
        end
    end

    % --- Final Assignment ---
    % Combine all indices and assign phase
    all_inds = vertcat(pore_linear_indices{:}); % [array] A 1D array vertically concatenating every single pore voxel index.
	
    % Overwrite the grain ID with the Pore ID. 
    % The Pore ID is deliberately set to N_grains + 1, so it is uniquely identifiable as 'not a grain'.	
    phase(all_inds) = N_grains + 1;

    % Statistics
    active_pores = sum(~cellfun(@isempty, pore_linear_indices)); % [integer] Counts how many pores actually intersected the domain and were inserted.
	
    nphase = N_grains + 1; % [integer] Total number of unique phases in the model (All grains + 1 global void phase).

    fprintf('Parallel Gas Void assignment: %d pores, %d voxels\n', ...
        N_pores, numel(all_inds));
    fprintf('Active pores (inside domain): %d out of %d (%.2f%%)\n', ...
        active_pores, N_pores, 100*active_pores/N_pores);
    fprintf('Average voxels per active pore: %.1f\n', ...
        numel(all_inds)/max(1,active_pores));
else
    nphase = N_grains;
    fprintf('No gas voids found to assign.\n');
end

fprintf('Phase assignment completed.\n');

%% 7. Assign Orientations -------------------------------------------------
% ==============================================================================
% Calculates and assigns 3D crystallographic orientations 
% (represented as Rodrigues vectors) to every single grain in the domain based 
% on the specified metallurgical texture type (Textured vs. Random).
% ==============================================================================

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
% ==============================================================================
% Converts Rodrigues vectors into MTEX orientation objects and calculates the 
% corresponding RGB color for each grain based on a standard IPF color key 
% relative to the sample's Z-axis.
% ==============================================================================

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
% ==============================================================================
% Exports the completed microstructure (grain/pore map) and orientation data to 
% disk. Supports both native MATLAB formats for internal review and universal 
% HDF5 formats for interfacing with external solid mechanics or FFT solvers.
% ==============================================================================

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
% ==============================================================================
% Generates a text file containing key statistics about the generated domain,
% including total grain/pore counts, physical dimensions, voxel resolution, 
% and calculated porosity percentages. Concludes by reporting total compute time.
% ==============================================================================

% Open file for writing
fid = fopen('output_summary.txt', 'w');

% Print summary of grains and pores found
fprintf(fid, 'Total number of Grains: %s, Pores: %s\n', ...
    regexprep(sprintf('%.0f',N_grains),'(\d)(?=(\d{3})+$)','$1,'), ...
    regexprep(sprintf('%.0f',nphase-N_grains),'(\d)(?=(\d{3})+$)','$1,'));

% Print physical size of the simulation domain and voxel resolution
fprintf(fid, 'Domain size: %.1f x %.1f x %.1f mm^3\n', ...
    cuboid_size(1), cuboid_size(2), cuboid_size(3));
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
% ==============================================================================
% Generates, formats, and saves a comprehensive suite of 2D and 3D figures to 
% visually and statistically validate the generated microstructure, porosity, 
% and crystallographic texture.
% ==============================================================================

if strcmp(micro_plot, 'on')

    %% --- 1. 2D Microplots (Cross-Sections) ---
    %==========================================================================
    % --- XY plane at Z = Nz (top surface) ---
    figure; clf
    s = squeeze(phase(:,:,Nz)); 
    imagesc(xvec, yvec, s');
    set(gca, 'YDir', 'normal'); axis equal tight
    colormap(cmap_ipf);
    set(gca, 'CLim', [1 size(cmap_ipf,1)]);
    xlabel('X [mm]'); ylabel('Y [mm]');
    title('Top Surface: XY Plane (Z = max)');
    set(gca, 'FontSize', 15, 'LineWidth', 1.1, 'Box', 'on')
    xticks(linspace(xmin, xmax, 6));
    yticks(linspace(ymin, ymax, 6));
    xlim([xmin xmax]);
    ylim([ymin ymax]);

    % --- YZ plane at X = 1 (left surface) ---
    figure; clf
    s = squeeze(phase(1,:,:)); 
    imagesc(yvec, zvec, s');
    set(gca, 'XDir', 'reverse');
    set(gca, 'YDir', 'normal');
    axis equal tight
    colormap(cmap_ipf);
    set(gca, 'CLim', [1 size(cmap_ipf,1)]);
    xlabel('Y [mm]'); ylabel('Z [mm]');
    title('Left Surface: YZ Plane (X = min)');
    set(gca, 'FontSize', 15, 'LineWidth', 1.1, 'Box', 'on')
    xticks(linspace(ymin, ymax, 6));
    yticks(linspace(zmin, zmax, 6));
    xlim([xmin xmax]);
    ylim([zmin zmax]);

    % --- XZ plane at Y = 1 (front surface) ---
    figure; clf
    s = squeeze(phase(:,1,:));
    imagesc(xvec, zvec, s');
    set(gca, 'YDir', 'normal'); axis equal tight
    colormap(cmap_ipf);
    set(gca, 'CLim', [1 size(cmap_ipf,1)]);
    xlabel('X [mm]'); ylabel('Z [mm]');
    title('Front Surface: XZ Plane (Y = min)');
    set(gca, 'FontSize', 15, 'LineWidth', 1.1, 'Box', 'on')
    xticks(linspace(xmin, xmax, 6));
    yticks(linspace(zmin, zmax, 6));
    xlim([xmin xmax]);
    ylim([zmin zmax]);

    %% --- 2. 3D Microstructure Scatter Plots ---

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
	
    %% --- 3. Statistical Histograms: 2D Section vs 3D Input ---

    % --- XY Plane Grain Size Histogram ---
    grain_diam_in_um = grain_diameters * 1e3;
    iz_mid = round(Nz/2);
    phase_xy = squeeze(phase(:, :, iz_mid));
    grain_ids = unique(phase_xy(:));
    grain_ids(grain_ids == N_grains+1) = [];
    grain_areas = zeros(numel(grain_ids),1);

    % Calculate equivalent circular diameter from voxel area in 2D
    for k = 1:numel(grain_ids)
        mask = (phase_xy == grain_ids(k));
        grain_areas(k) = sum(mask(:)) * dx * dx;
    end
    grain_diam_xy = 2*sqrt(grain_areas/pi);
    grain_diam_xy_um = grain_diam_xy * 1e3;

    % Create figure for histogram
    figure; clf
    data = grain_diam_xy_um;
	
    % Automatic binning using Doane's Formula (better for highly skewed lognormal data)	
    skewness_val = skewness(data);
    num_bin = ceil(1 + log2(numel(data)) + log2(1 + abs(skewness_val) / sqrt(6/numel(data))));
    h1 = histogram(grain_diam_xy_um, num_bin, 'Normalization','pdf', ...
        'FaceAlpha',0.3, 'FaceColor',[0.2 0.8 0.2],'EdgeColor','k', ...
        'DisplayName','Grain size (XY Plane)');
    hold on;

    % Fit Lognormal Distribution to 2D section data
    pd_xy = fitdist(grain_diam_xy_um, 'Lognormal');
    xfit = linspace(0, max(grain_diam_xy_um)*1.05, 200);
    yfit2 = pdf(pd_xy, xfit);
    plot(xfit, yfit2, 'b-', 'LineWidth', 3, 'DisplayName', 'Lognormal fit');
    hold off

    xl = prctile(grain_diam_xy_um, [0 100]);
    xlim([min(0,xl(1)), xl(2)*1.1]);
    ylim([0, max(h1.Values)*1.2]);
    xlabel('Grain diameter [\mum]', 'FontSize',15);
    ylabel('Probability density', 'FontSize',15);
    title('Grain Size: Input Distribution in XY Section (Z = max)', 'FontSize',15);
    set(gca, 'FontSize',15); legend('show'); hold off
	
    % Annotate with statistics
    mean_fit = exp(pd_xy.mu + 0.5 * pd_xy.sigma^2);
    std_fit  = sqrt((exp(pd_xy.sigma^2) - 1) * exp(2*pd_xy.mu + pd_xy.sigma^2));
    xl_text = xlim; yl_text = ylim;
    x_text = xl_text(2) - 0.35 * (xl_text(2) - xl_text(1));
    y_text = yl_text(1) + 0.5 * (yl_text(2) - yl_text(1));

    text(x_text, y_text, { ...
        sprintf('Grain size: mean = %.2f \\mum', mean_fit), ...
        sprintf('Grain size: std = %.2f \\mum', std_fit), ...
        sprintf('lognormal: \\mu = %.2f, \\sigma = %.2f', pd_xy.mu, pd_xy.sigma)}, ...
        'BackgroundColor','w', 'FontSize',15, 'HorizontalAlignment','left');


    % --- 3D Input Grain Size Histogram ---
    figure; clf; hold on;
    grain_diam_in_um = grain_diameters * 1e3;	
    data_3d = grain_diam_in_um;
    skewness_val = skewness(data_3d);
    num_bin = ceil(1 + log2(numel(data_3d)) + log2(1 + abs(skewness_val) / sqrt(6/numel(data_3d))));

    h2 = histogram(grain_diam_in_um, num_bin, 'Normalization','pdf', ...
        'FaceAlpha',0.6, 'FaceColor',[0.2 0.8 0.2],'EdgeColor','k', ...
        'DisplayName','3D input');

    pd_in = fitdist(grain_diam_in_um, 'Lognormal');
    xfit = linspace(0, max(grain_diam_in_um), 200);
    yfit = pdf(pd_in, xfit);
    plot(xfit, yfit, 'r-', 'LineWidth', 2, 'DisplayName', '3D Lognormal fit');
    hold off

    xl = prctile(grain_diam_in_um, [0 100]);
    xlim([max(0,xl(1)), xl(2)*1.05]);
    yl = ylim;
    ylim([0, max(h2.Values)*1.2]);
    xlabel('Grain diameter [\mum]', 'FontSize',15);
    ylabel('Probability density', 'FontSize',15);
    title('Grain Size: 3D Input Distribution', 'FontSize',15);
    set(gca, 'FontSize',15); legend('show'); hold off

    x_text = xl(2) - 0.35 * (xl(2) - xl(1));  
    y_text = yl(1) + 0.5 * (yl(2) - yl(1));   
    mean_fit = exp(pd_in.mu + 0.5 * pd_in.sigma^2);
    std_fit  = sqrt((exp(pd_in.sigma^2) - 1) * exp(2*pd_in.mu + pd_in.sigma^2));

    text(x_text, y_text, { ...
        sprintf('Grain size: mean = %.2f \\mum', mean_fit), ...
        sprintf('Grain size: std = %.2f \\mum', std_fit), ...
        sprintf('lognormal: \\mu = %.2f, \\sigma = %.2f', pd_in.mu, pd_in.sigma)}, ...
        'BackgroundColor', 'w', 'FontSize', 15);


    %% --- 4. Aspect Ratio Histograms ---
    % -- XZ PLANE Aspect Ratio (2D) --
    iy_mid = round(Ny/2);
    phase_xz = squeeze(phase(:, iy_mid, :));
    grain_ids_xz = unique(phase_xz(:));
    grain_ids_xz(grain_ids_xz == N_grains+1) = [];
    ratios_xz = zeros(numel(grain_ids_xz),1);

    for k = 1:numel(grain_ids_xz)
        id = grain_ids_xz(k);
        [ix, iz] = find(phase_xz == id);
        if isempty(ix)
            ratios_xz(k) = nan;
        else
            x_extent = max(xvec(ix)) - min(xvec(ix));
            z_extent = max(zvec(iz)) - min(zvec(iz));
            if x_extent==0 || z_extent==0
                ratios_xz(k) = nan;
            else
                ratios_xz(k) = z_extent / x_extent;
            end
        end
    end
    ratios_xz = ratios_xz(~isnan(ratios_xz) & ratios_xz > 0);
    input_aspect_xz = aspect_ratios;
	
    data = ratios_xz;
    skewness_val = skewness(data);
    num_bin = ceil(1 + log2(numel(data)) + log2(1 + abs(skewness_val) / sqrt(6/numel(data))));

    figure; clf
    h1 = histogram(ratios_xz, num_bin, 'Normalization','pdf', ...
        'FaceAlpha',0.3,'FaceColor',[0.9 0.4 0.2],'EdgeColor','k', ...
        'DisplayName','Aspect ratio (XZ Plane)');
    hold on;

    xfit = linspace(min([ratios_xz; input_aspect_xz]), ...
        max([ratios_xz; input_aspect_xz]), 500);
    pd_2d = fitdist(ratios_xz, 'Lognormal');
    yfit2 = pdf(pd_2d, xfit);
    plot(xfit, yfit2, 'b-', 'LineWidth',3,'DisplayName','Lognormal fit');
    hold off;

    xl = prctile(data, [0 100]);
    xlim([max(0, xl(1)), xl(2)*1.05]);
    yl = ylim;
    ylim([0, max(h1.Values)*1.2]);
    xlabel('Length/Width ratio (XZ plane)', 'FontSize',15);
    ylabel('Probability density', 'FontSize',15);
    title('Aspect ratio: XZ plane (Y = min)', 'FontSize',15);
    set(gca, 'FontSize',15); legend('show'); hold off

    mean_fit = exp(pd_2d.mu + 0.5 * pd_2d.sigma^2);
    std_fit  = sqrt((exp(pd_2d.sigma^2) - 1) * exp(2*pd_2d.mu + pd_2d.sigma^2));
    x_text = xl(2) - 0.35 * (xl(2) - xl(1));  
    y_text = yl(1) + 0.5 * (yl(2) - yl(1));

    text(x_text, y_text, { ...
        sprintf('Aspect ratio: mean = %.2f', mean_fit), ...
        sprintf('Aspect ratio: std = %.2f', std_fit), ...
        sprintf('lognormal: \\mu = %.2f, \\sigma = %.2f', pd_2d.mu, pd_2d.sigma)}, ...
        'BackgroundColor', 'w', 'FontSize', 15);

    % -- 3D Input Aspect Ratio Histogram (XZ) --
    figure; clf
    input_aspect_xz = aspect_ratios;
    skewness_val = skewness(input_aspect_xz);
    num_bin = ceil(1 + log2(numel(input_aspect_xz)) + log2(1 + abs(skewness_val) / sqrt(6/numel(input_aspect_xz))));

    hold on;
    h2 = histogram(input_aspect_xz, num_bin, 'Normalization','pdf', ...
        'FaceAlpha',0.6,'FaceColor',[0.9 0.4 0.2],'EdgeColor','k', ...
        'DisplayName','3D input');
    pd_in = fitdist(input_aspect_xz, 'Lognormal');
    xfit = linspace(min(input_aspect_xz), max(input_aspect_xz), 500);
    yfit = pdf(pd_in, xfit);
    plot(xfit, yfit, 'r-', 'LineWidth',2, 'DisplayName','3D Lognormal fit');
    hold off;

    xl = prctile(input_aspect_xz, [0 100]);
    xlim([max(0,xl(1)), xl(2)*1.1]);
    yl = ylim;
    ylim([0, max([h1.Values, yfit])*1.2]);
    xlabel('Length/Width ratio (XZ plane)', 'FontSize',15);
    ylabel('Probability density', 'FontSize',15);
    title('Aspect ratio: 3D Input', 'FontSize',15);
    set(gca, 'FontSize',15); legend('show'); hold off

    mean_fit = exp(pd_in.mu + 0.5 * pd_in.sigma^2);
    std_fit  = sqrt((exp(pd_in.sigma^2) - 1) * exp(2*pd_in.mu + pd_in.sigma^2));
    x_text = xl(2) - 0.25 * (xl(2) - xl(1));  
    y_text = yl(1) + 0.5 * (yl(2) - yl(1)); 

    text(x_text, y_text, { ...
        sprintf('Aspect ratio: mean = %.2f', mean_fit), ...
        sprintf('Aspect ratio: std = %.2f', std_fit), ...
        sprintf('lognormal: \\mu = %.2f, \\sigma = %.2f', pd_in.mu, pd_in.sigma)}, ...
        'BackgroundColor', 'w', 'FontSize', 15);


    % -- YZ PLANE Aspect Ratio Analysis --
    ix_mid = round(Nx/2);
    phase_yz = squeeze(phase(ix_mid, :, :));
    grain_ids_yz = unique(phase_yz(:));
    grain_ids_yz(grain_ids_yz == N_grains+1) = [];
    ratios_yz = zeros(numel(grain_ids_yz),1);
    for k = 1:numel(grain_ids_yz)
        id = grain_ids_yz(k);
        [iy, iz] = find(phase_yz == id);
        if isempty(iy)
            ratios_yz(k) = nan;
        else
            y_extent = max(yvec(iy)) - min(yvec(iy));
            z_extent = max(zvec(iz)) - min(zvec(iz));
            if y_extent==0 || z_extent==0
                ratios_yz(k) = nan;
            else
                ratios_yz(k) = z_extent / y_extent;
            end
        end
    end
    ratios_yz = ratios_yz(~isnan(ratios_yz) & ratios_yz > 0);
    input_aspect_yz = aspect_ratios;

    figure; clf
    skewness_val = skewness(ratios_yz);
    num_bin = ceil(1 + log2(numel(ratios_yz)) + log2(1 + abs(skewness_val) / sqrt(6/numel(ratios_yz))));
    h1 = histogram(ratios_yz, num_bin, 'Normalization','pdf', ...
        'FaceAlpha',0.3,'FaceColor',[0.2 0.4 0.9],'EdgeColor','k', ...
        'DisplayName','Aspect ratio (YZ Plane)');
    hold on;

    xfit = linspace(min([ratios_yz; input_aspect_yz]), ...
        max([ratios_yz; input_aspect_yz]), 500);
    pd_2d = fitdist(ratios_yz, 'Lognormal');
    yfit2 = pdf(pd_2d, xfit);
    plot(xfit, yfit2, 'b-', 'LineWidth',3,'DisplayName','Lognormal fit');

    xl = prctile(ratios_yz, [0 100]);
    xlim([max(0,xl(1)), xl(2)*1.1]);
    yl = ylim;
    ylim([0, max([h1.Values, yfit2])*1.2]);

    xlabel('Length/Width ratio (YZ plane)', 'FontSize',15);
    ylabel('Probability density', 'FontSize',15);
    title('Aspect ratio: YZ plane (X = min)', 'FontSize',15);
    set(gca, 'FontSize',15); legend('show'); hold off
	
    mean_fit = exp(pd_2d.mu + 0.5 * pd_2d.sigma^2);
    std_fit  = sqrt((exp(pd_2d.sigma^2) - 1) * exp(2*pd_2d.mu + pd_2d.sigma^2));
    x_text = xl(2) - 0.35 * (xl(2) - xl(1)); 
    y_text = yl(1) + 0.5 * (yl(2) - yl(1));

    text(x_text, y_text, { ...
        sprintf('Aspect ratio: mean = %.2f', mean_fit), ...
        sprintf('Aspect ratio: std = %.2f', std_fit), ...
        sprintf('lognormal: \\mu = %.2f, \\sigma = %.2f', pd_2d.mu, pd_2d.sigma)}, ...
        'BackgroundColor','w','FontSize',15);


    % -- 3D Input Aspect Ratio Histogram (YZ) --
    figure; clf
    input_aspect_yz = aspect_ratios;
    skewness_val = skewness(input_aspect_yz);
    num_bin = ceil(1 + log2(numel(input_aspect_yz)) + log2(1 + abs(skewness_val) / sqrt(6/numel(input_aspect_yz))));

    hold on;
    h2 = histogram(input_aspect_yz, num_bin, 'Normalization','pdf', ...
        'FaceAlpha',0.6,'FaceColor',[0.2 0.4 0.9],'EdgeColor','k', ...
        'DisplayName','3D input');
    pd_in = fitdist(input_aspect_yz, 'Lognormal');
    xfit = linspace(min([ratios_yz; input_aspect_yz]), ...
        max([ratios_yz; input_aspect_yz]), 500);
    yfit = pdf(pd_in, xfit);
    plot(xfit, yfit, 'r-', 'LineWidth',2, 'DisplayName','3D Lognormal fit');
    hold off;

    xl = prctile(input_aspect_yz, [0 100]);
    xlim([max(0,xl(1)), xl(2)*1.05]);
    yl = ylim;
    ylim([0, max([h2.Values, yfit])*1.2]);
    xlabel('Length/Width ratio (YZ plane)', 'FontSize',15);
    ylabel('Probability density', 'FontSize',15);
    title('Aspect ratio: 3D Input', 'FontSize',15);
    set(gca, 'FontSize',15); legend('show'); hold off

    mean_fit = exp(pd_in.mu + 0.5 * pd_in.sigma^2);
    std_fit  = sqrt((exp(pd_in.sigma^2) - 1) * exp(2*pd_in.mu + pd_in.sigma^2));
    x_text = xl(2) - 0.3 * (xl(2) - xl(1));  
    y_text = yl(1) + 0.5 * (yl(2) - yl(1));

    text(x_text, y_text, { ...
        sprintf('Aspect ratio: mean = %.2f', mean_fit), ...
        sprintf('Aspect ratio: std = %.2f', std_fit), ...
        sprintf('lognormal: \\mu = %.2f, \\sigma = %.2f', pd_in.mu, pd_in.sigma)}, ...
        'BackgroundColor','w','FontSize',15);


    %% --- 5. Advanced Isosurface Rendering ---
    % Renders actual voxel shapes (smoothed polygons) rather than scatter points
    pore_color     = [0 0 1];
    pore_id        = N_grains + 1;

    figure; clf; hold on;

    % Plot a subset of Grains
    n_show = min(num_grains_iso, N_grains);
    cmap_grains = lines(max(64, n_show));
    if n_show > size(cmap_grains,1)
        cmap_grains = parula(n_show);
    end

    fprintf('Plotting %d largest grains...\n', n_show);

    for k = 1:n_show
        grain_mask = (phase == k);

        if any(grain_mask(:))
            fv = isosurface(xvec, yvec, zvec, smooth3(double(grain_mask),'box',3), 0.5);

            if ~isempty(fv.vertices)
                fv_red = reducepatch(fv, 0.15);
                patch(fv_red, 'FaceColor', cmap_grains(k,:), ...
                    'EdgeColor', 'none', 'FaceAlpha', 0.38, ...
                    'SpecularStrength', 0.2);
            end
        end
    end

    % Plot Selective Pores
    fprintf('Analyzing pore geometry for plotting...\n');
    pore_mask_all = (phase == pore_id);

    if any(pore_mask_all(:))
        CC = bwconncomp(pore_mask_all, 6);
        stats = regionprops(CC, 'Area', 'PixelIdxList');
        [~, sort_idx] = sort([stats.Area], 'descend');
        n_plot_actual = min(max_pores_show, length(sort_idx));

        selective_mask = false(size(pore_mask_all));
        for i = 1:n_plot_actual
            idx = sort_idx(i);
            selective_mask(stats(idx).PixelIdxList) = true;
        end

        fv_pores = isosurface(xvec, yvec, zvec, double(selective_mask), 0.5);

        if ~isempty(fv_pores.vertices)
            patch(fv_pores, ...
                'FaceColor', pore_color, ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 1.0, ...
                'DiffuseStrength', 0.8, ...
                'SpecularStrength', 0.5, ...
                'SpecularExponent', 10);
        end
    else
        n_plot_actual = 0;
    end


    view([-42.71 22.97]);
    axis equal tight; grid on; box on
    set(gca, 'Color', [0.9 0.9 0.9]);
    camlight('headlight');
    lighting gouraud;

    hx = xlabel('X_1 (mm)'); set(hx, 'Rotation', 21);
    hy = ylabel('X_2 (mm)'); set(hy, 'Rotation', -21);
    zlabel('X_3 (mm)');
    title({sprintf('Microstructure: Top %d Pores + %d Grains', n_plot_actual, n_show)}, ...
        'FontSize', 15);

    xticks(linspace(xmin, xmax, 5));
    yticks(linspace(ymin, ymax, 5));
    zticks(linspace(zmin, zmax, 5));
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    zlim([zmin,zmax]);

    set(gca, 'FontSize', 15, 'LineWidth', 1.3);
    hAxes = findobj(gcf,"Type","axes");
    if ~isempty(hAxes)
        hAxes.BoxStyle = "full";
        hAxes.ClippingStyle = "3dbox";
    end
    hold off;


    %% --- 6. Crystallographic Pole Figures (MTEX) ---

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

    % --- Plot Inverse Pole Figure (IPF) ---
    figure; clf;
	
    % Project standard sample directions (X,Y,Z) back onto the fundamental crystal triangle
    x_dir = xvector; 
    y_dir = yvector; 
    z_dir = zvector;
    plotIPDF(odf, [x_dir, y_dir, z_dir], 'antipodal');


elseif strcmp(micro_plot, 'off')
    disp('micro_plot_disabled')
end

%% ---- END ----
