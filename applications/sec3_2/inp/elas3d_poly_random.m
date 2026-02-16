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
%   - MTEX Toolbox (Required for orientation assignment and pole figures)
%   - Parallel Computing Toolbox (Required for spatial filtering acceleration)
%
% Inputs:
%   - Configuration parameters (Domain size, Grid resolution, Grain statistics)
%   - 'Defect_data.xlsx': Experimental defect data containing pore centers and radii
%
% Outputs:
%   - 'input_structure_poly.h5': HDF5 file containing the 3D phase map and orientations
%   - 'xct_poly_params_SS316L.mat': Saved workspace parameters
%   - 'output_summary.txt': Text file containing volume fractions and mesh statistics
%   - Visualization figures: 3D scatter surfaces, smooth isosurfaces, and pole figures
%
% Core Workflow:
%   1. Define computational domain and voxel discretization.
%   2. Generate grain seeds using lognormal size and aspect ratio distributions.
%   3. Assign grain phases via highly optimized, parallelized spatial filtering.
%   4. Map explicit defect/pore regions into the voxel grid.
%   5. Assign crystallographic orientations (Uniform Random or Fiber texture).
%   6. Export binary HDF5 structures and render diagnostics.
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
% ------------------------------------------------------------------------------
format compact; clear; close all;
tStart = tic;
%% 1. Configuration & Parameters -------------------------------------------------------------
max_len = 0.2;      
cube_len = 0.2;    
grid_num = 800;     
dx = cube_len/grid_num;
xmin = -0.1;
ymin = -0.1;
zmin = -0.1;

% --- Grain Statistics ---
mean_d = 0.005;        
std_d  = 0.10*mean_d;  

% Aspect ratio (length/width)
mean_aspect = 1.000;            
std_aspect = 0.1 * mean_aspect;

% --- Output Settings Variables ---
phase_data_type = 'uint32';

% Selctive plot of grains
stride = 1;
microstructure3d_type = 'surface';
max_pores_show = 1;
num_grains_iso = 10;
K_grain = 2000; 

% --- Pore Settings ---
pore = 'on';
input_pore_file = 'Defect_data.xlsx';
%% 2. Pore Data -----------------------------------------------------------
if strcmp(pore, 'on')
    data = readmatrix(input_pore_file);
    x_pore       = data(:,9);     
    y_pore       = data(:,10);     
    z_pore       = data(:,11);     
    d = data(:,8);     
    r = d/2;          
    V_pore = data(:,7);    
    D_Feret_max = data(:,12);
    N_pores = numel(x_pore);

elseif strcmp(pore, 'off')
    x_pore = [];
    y_pore = [];
    z_pore = [];
    d = [];
    r = [];
    V_pore = [];
    D_Feret_max = [];
    N_pores = 0;
    warning('No pore file found. Proceeding with no pores.');
    % exit;
end

%% 3. Computational Domain ------------------------------------------------
cube_len = min(cube_len, max_len); 

xmax = xmin + cube_len;
ymax = ymin + cube_len;
zmax = zmin + cube_len;

cuboid_size = [xmax-xmin, ymax-ymin, zmax-zmin];

xvec = xmin:dx:xmax;   
Nx = length(xvec);

yvec = ymin:dx:ymax;   
Ny = length(yvec);

zvec = zmin:dx:zmax;   
Nz = length(zvec);

[X, Y, Z] = ndgrid(xvec, yvec, zvec);

domain_vol = prod(cuboid_size);


%% 4. Lognormal Grain Size and Aspect Ratio --------------------------------
mean_grain_vol = pi*(mean_d/2)^2 * (mean_aspect*mean_d);
N_grains = ceil(domain_vol / mean_grain_vol);

sigma_log = sqrt(log((std_d/mean_d)^2 + 1));
mu_log = log(mean_d) - 0.5*sigma_log^2;
grain_diameters = lognrnd(mu_log, sigma_log, N_grains, 1); 
r_xy = grain_diameters / 2;

sigma_log_a = sqrt(log((std_aspect/mean_aspect)^2 + 1));
mu_log_a    = log(mean_aspect) - 0.5*sigma_log_a^2;
aspect_ratios = lognrnd(mu_log_a, sigma_log_a, N_grains, 1); 
r_z  = aspect_ratios .* r_xy;

grain_seeds = [rand(N_grains,1)*cuboid_size(1)+xmin, ...
    rand(N_grains,1)*cuboid_size(2)+ymin, ...
    rand(N_grains,1)*cuboid_size(3)+zmin];

%% 5.Anisotropic Assignment -------------------------------------
% 5-1: Print header for upcoming block assignment step
fprintf('\n=== PARALLEL + SPATIAL FILTERING MICROSTRUCTURE GENERATION ===\n');
main_timer = tic;

if isempty(gcp('nocreate'))
    fprintf('Starting parallel pool...\n');
    pool_timer = tic;
    num_workers = feature('numcores');
    parpool('local', num_workers);
    fprintf('Parallel pool started in %.1fs with %d workers\n', ...
        toc(pool_timer), gcp().NumWorkers);
else
    fprintf('Using existing parallel pool with %d workers\n', gcp().NumWorkers);
end

phase = zeros(Nx, Ny, Nz, phase_data_type);
inv_r_xy2 = single(1 ./ (r_xy.^2));
inv_r_z2 = single(1 ./ (r_z.^2));
grain_seeds_single = single(grain_seeds);

% 5-2 Build Spatial Acceleration Structure
fprintf('\nBuilding spatial acceleration structure...\n');
accel_timer = tic;

grid_res = 15;
x_edges = linspace(min(xvec), max(xvec), grid_res+1);
y_edges = linspace(min(yvec), max(yvec), grid_res+1);
z_edges = linspace(min(zvec), max(zvec), grid_res+1);

max_influence = max([r_xy(:); r_z(:)]) * 1.5;
total_cells = grid_res^3;

spatial_grid_temp = cell(total_cells, 1);

parfor linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);

    cell_xmin = x_edges(i); cell_xmax = x_edges(i+1);
    cell_ymin = y_edges(j); cell_ymax = y_edges(j+1);
    cell_zmin = z_edges(k); cell_zmax = z_edges(k+1);

    relevant_grains = find(...
        grain_seeds_single(:,1) >= (cell_xmin - max_influence) & ...
        grain_seeds_single(:,1) <= (cell_xmax + max_influence) & ...
        grain_seeds_single(:,2) >= (cell_ymin - max_influence) & ...
        grain_seeds_single(:,2) <= (cell_ymax + max_influence) & ...
        grain_seeds_single(:,3) >= (cell_zmin - max_influence) & ...
        grain_seeds_single(:,3) <= (cell_zmax + max_influence));

    spatial_grid_temp{linear_idx} = relevant_grains;
end

spatial_grid = cell(grid_res, grid_res, grid_res);
for linear_idx = 1:total_cells
    [i, j, k] = ind2sub([grid_res, grid_res, grid_res], linear_idx);
    spatial_grid{i,j,k} = spatial_grid_temp{linear_idx};
end
clear spatial_grid_temp;

fprintf('Spatial acceleration structure built in %.1fs\n', toc(accel_timer));

% 5-3 Setup Parallel Block Processing
block_size = 50;
n_blocks_x = ceil(Nx / block_size);
n_blocks_y = ceil(Ny / block_size);
n_blocks_z = ceil(Nz / block_size);
total_blocks = n_blocks_x * n_blocks_y * n_blocks_z;

fprintf('Processing %d blocks (%dx%dx%d) in parallel:\n', total_blocks, n_blocks_x, n_blocks_y, n_blocks_z);

num_workers = gcp().NumWorkers;
blocks_per_batch = max(1, floor(total_blocks / (num_workers * 4)));
n_batches = ceil(total_blocks / blocks_per_batch);

fprintf('Using %d batches of %d blocks each for progress tracking\n', n_batches, blocks_per_batch);

completed_blocks = 0;
total_grains_processed = 0;
blocks_with_grains = 0;
batch_times = [];

% 5-4 Process Batches with Progress Tracking
fprintf('\nStarting parallel processing with real-time progress:\n');
process_timer = tic;

for batch_idx = 1:n_batches
    batch_timer = tic;

    batch_start = (batch_idx - 1) * blocks_per_batch + 1;
    batch_end = min(batch_idx * blocks_per_batch, total_blocks);
    current_batch_size = batch_end - batch_start + 1;

    block_list_batch = zeros(current_batch_size, 4);
    for local_idx = 1:current_batch_size
        block_num = batch_start + local_idx - 1;
        bz = ceil(block_num / (n_blocks_x * n_blocks_y));
        remaining = block_num - (bz-1) * n_blocks_x * n_blocks_y;
        by = ceil(remaining / n_blocks_x);
        bx = remaining - (by-1) * n_blocks_x;
        block_list_batch(local_idx, :) = [bx, by, bz, block_num];
    end

    batch_results = cell(current_batch_size, 1);
    batch_indices = cell(current_batch_size, 1);
    batch_grain_counts = zeros(current_batch_size, 1);

    xvec_local = xvec;
    yvec_local = yvec;
    zvec_local = zvec;
    x_edges_local = x_edges;
    y_edges_local = y_edges;
    z_edges_local = z_edges;


    parfor local_idx = 1:current_batch_size
        bx = block_list_batch(local_idx, 1);
        by = block_list_batch(local_idx, 2);
        bz = block_list_batch(local_idx, 3);

        ix1 = (bx-1)*block_size + 1; ix2 = min(bx*block_size, Nx);
        iy1 = (by-1)*block_size + 1; iy2 = min(by*block_size, Ny);
        iz1 = (bz-1)*block_size + 1; iz2 = min(bz*block_size, Nz);

        block_xmin = xvec_local(ix1); block_xmax = xvec_local(ix2);
        block_ymin = yvec_local(iy1); block_ymax = yvec_local(iy2);
        block_zmin = zvec_local(iz1); block_zmax = zvec_local(iz2);

        relevant_grains = find(...
            grain_seeds_single(:,1) >= (block_xmin - max_influence) & ...
            grain_seeds_single(:,1) <= (block_xmax + max_influence) & ...
            grain_seeds_single(:,2) >= (block_ymin - max_influence) & ...
            grain_seeds_single(:,2) <= (block_ymax + max_influence) & ...
            grain_seeds_single(:,3) >= (block_zmin - max_influence) & ...
            grain_seeds_single(:,3) <= (block_zmax + max_influence));

        if ~isempty(relevant_grains)
            [Xb, Yb, Zb] = ndgrid(single(xvec_local(ix1:ix2)), single(yvec_local(iy1:iy2)), single(zvec_local(iz1:iz2)));

            block_phase = zeros(size(Xb), 'uint32');
            min_dist2 = inf(size(Xb), 'single');

            for k = relevant_grains'
                dX = Xb - grain_seeds_single(k,1);
                dY = Yb - grain_seeds_single(k,2);
                dZ = Zb - grain_seeds_single(k,3);

                current_dist2 = (dX.^2) * inv_r_xy2(k) + (dY.^2) * inv_r_xy2(k) + (dZ.^2) * inv_r_z2(k);

                update_mask = current_dist2 < min_dist2;
                block_phase(update_mask) = k;
                min_dist2(update_mask) = current_dist2(update_mask);
            end

            batch_results{local_idx} = block_phase;
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = length(relevant_grains);
        else
            batch_results{local_idx} = [];
            batch_indices{local_idx} = [ix1, ix2, iy1, iy2, iz1, iz2];
            batch_grain_counts(local_idx) = 0;
        end
    end

    batch_grains_processed = 0;
    batch_blocks_with_grains = 0;

    for local_idx = 1:current_batch_size
        if ~isempty(batch_results{local_idx})
            indices = batch_indices{local_idx};
            ix1 = indices(1); ix2 = indices(2);
            iy1 = indices(3); iy2 = indices(4);
            iz1 = indices(5); iz2 = indices(6);

            phase(ix1:ix2, iy1:iy2, iz1:iz2) = cast(batch_results{local_idx}, phase_data_type);
            batch_blocks_with_grains = batch_blocks_with_grains + 1;
        end
        batch_grains_processed = batch_grains_processed + batch_grain_counts(local_idx);
    end

    completed_blocks = completed_blocks + current_batch_size;
    total_grains_processed = total_grains_processed + batch_grains_processed;
    blocks_with_grains = blocks_with_grains + batch_blocks_with_grains;

    batch_time = toc(batch_timer);
    batch_times(end+1) = batch_time;

    progress_pct = 100 * completed_blocks / total_blocks;
    elapsed_total = toc(main_timer);
    avg_grains_per_block = batch_grains_processed / current_batch_size;

    if batch_idx > 1
        avg_batch_time = mean(batch_times);
        remaining_batches = n_batches - batch_idx;
        eta_seconds = remaining_batches * avg_batch_time;
    else
        eta_seconds = batch_time * (n_batches - 1);
    end

    fprintf('Batch %d/%d: Blocks %d-%d (%.1f%%) | %.1f grains/block avg | Elapsed: %.1fs | ETA: %.1fs (%.1f min)\n', ...
        batch_idx, n_batches, batch_start, batch_end, progress_pct, avg_grains_per_block, elapsed_total, eta_seconds, eta_seconds/60);

    if mod(batch_idx, max(1, floor(n_batches/5))) == 0 || batch_idx == n_batches
        overall_avg_grains = total_grains_processed / completed_blocks;
        processing_rate = completed_blocks / elapsed_total;
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

unassigned = sum(phase(:) == 0);
if unassigned > 0
    fprintf('\nWARNING: %d unassigned voxels detected!\n', unassigned);
else
    fprintf('\nSUCCESS: All %s voxels assigned successfully!\n', regexprep(sprintf('%.0f',double(Nx*Ny*Nz)),'(\d)(?=(\d{3})+$)','$1,'));
end

fprintf('\nParallel + Spatial filtering microstructure generation completed!\n');



%% 6. Overwrite With Pore Regions -----------------------------------------
sz_phase = size(phase);
Nx = sz_phase(1); Ny = sz_phase(2); Nz = sz_phase(3);

pore_linear_indices = cell(N_pores, 1);

if ~isempty(x_pore)
    fprintf('Inserting %d Gas Voids (Spherical)... \n', N_pores);

    parfor k = 1:N_pores

        r_pore = d(k) / 2.0;

        if r_pore <= 0 || isnan(r_pore)
            continue;
        end

        buffer = 1.1;
        max_r = r_pore * buffer;

        x_min = x_pore(k) - max_r;
        x_max = x_pore(k) + max_r;
        y_min = y_pore(k) - max_r;
        y_max = y_pore(k) + max_r;
        z_min = z_pore(k) - max_r;
        z_max = z_pore(k) + max_r;

        idx_x = find(xvec >= x_min & xvec <= x_max);
        idx_y = find(yvec >= y_min & yvec <= y_max);
        idx_z = find(zvec >= z_min & zvec <= z_max);

        if isempty(idx_x) || isempty(idx_y) || isempty(idx_z)
            continue;
        end

        [X_sub, Y_sub, Z_sub] = ndgrid(xvec(idx_x), yvec(idx_y), zvec(idx_z));

        Xc = X_sub - x_pore(k);
        Yc = Y_sub - y_pore(k);
        Zc = Z_sub - z_pore(k);

        mask_sub = (Xc.^2 + Yc.^2 + Zc.^2) <= r_pore^2;

        [ii, jj, kk] = ind2sub(size(mask_sub), find(mask_sub));

        if ~isempty(ii)
            global_x = idx_x(ii);
            global_y = idx_y(jj);
            global_z = idx_z(kk);

            pore_linear_indices{k} = sub2ind(sz_phase, ...
                double(global_x(:)), ...
                double(global_y(:)), ...
                double(global_z(:)));
        end
    end

    all_inds = vertcat(pore_linear_indices{:});
    phase(all_inds) = N_grains + 1;

    active_pores = sum(~cellfun(@isempty, pore_linear_indices));
    nphase = N_grains + 1;

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
ori_vec = zeros(3*N_grains, 1);

fprintf('Generating Random Texture \n');
for k = 1:N_grains
    u1 = rand(); u2 = rand(); u3 = rand();

    q_w = sqrt(1-u1) * sin(2*pi*u2);
    q_x = sqrt(1-u1) * cos(2*pi*u2);
    q_y = sqrt(u1)   * sin(2*pi*u3);
    q_z = sqrt(u1)   * cos(2*pi*u3);

    if abs(q_w) < 1e-8
        rod_vec = [q_x q_y q_z] / 1e-9;
    else
        rod_vec = [q_x q_y q_z] / q_w;
    end

    idx = 3*(k-1)+1;
    ori_vec(idx:idx+2) = rod_vec;
end

%% 8. --- Preparation for IPF coloring ---
cs = crystalSymmetry('m-3m');
rods = reshape(ori_vec(1:3*N_grains), [3, N_grains])';
ori_mtex = orientation('rodrigues', rods, cs);
ipf_key = ipfColorKey(cs);
ipf_key.inversePoleFigureDirection = vector3d.Z;
cmap_ipf = ipf_key.orientation2color(ori_mtex);
	
%% 9. Save To HDF5 --------------------------------------------------------
n_phases = N_grains + 1;
save('xct_poly_params_SS316L.mat','xvec','yvec','zvec','Nx','Ny','Nz', ...
    'cube_len','N_grains','n_phases');
save('xct_poly_params_SS316L_full.mat','-v7.3')

h5filename = 'input_structure_poly.h5';
if isfile(h5filename)
    delete(h5filename);
    fprintf('Previous file "%s" deleted.\n', h5filename);
end

pix_flat = int32(phase(:));
h5create(h5filename, '/pix', [length(pix_flat) 1], 'Datatype', 'int32');
h5write(h5filename, '/pix', pix_flat);
h5writeatt(h5filename, '/pix', 'dimensions', int32([Nx Ny Nz]));
h5create(h5filename, '/orientation', size(ori_vec), 'Datatype', 'double');
h5write(h5filename, '/orientation', ori_vec);
h5writeatt(h5filename, '/orientation', 'dimensions', int32(size(ori_vec)));

disp(['Saved HDF5 input structure to ', h5filename]);


%% 10. Output Summary
fid = fopen('output_summary.txt', 'w');

fprintf(fid, 'Total number of Grains: %s, Pores: %s\n', ...
    regexprep(sprintf('%.0f',N_grains),'(\d)(?=(\d{3})+$)','$1,'), ...
    regexprep(sprintf('%.0f',nphase-N_grains),'(\d)(?=(\d{3})+$)','$1,'));
	
fprintf(fid, 'Domain: %d x %d x %d = %s voxels\n', ...
    Nx, Ny, Nz, regexprep(sprintf('%.0f',Nx*Ny*Nz),'(\d)(?=(\d{3})+$)','$1,'));
fprintf('Grains: %s\n', regexprep(sprintf('%.0f',N_grains),'(\d)(?=(\d{3})+$)','$1,'));

num_vox_grain = sum( phase(:) >= 1 & phase(:) <= N_grains );
num_vox_pore = sum( phase(:) == (N_grains+1) );
num_vox_total = numel(phase);

fprintf(fid, 'Grain voxels: %d (%.3f%%)\n', num_vox_grain, 100*num_vox_grain/num_vox_total);
fprintf(fid, 'Pore  voxels: %d (%.3f%%)\n', num_vox_pore,  100*num_vox_pore/num_vox_total);
fprintf(fid, 'Total voxels: %d\n', num_vox_total);


fprintf(fid,'Mean grain size (XY): %.3f um\n',mean(grain_diameters)*1e3);
fprintf(fid,'Mean aspect ratio (Z/XY): %.3f\n',mean(aspect_ratios));
fprintf(fid,'Done! Phase voxel mesh saved for full cuboid domain with defects.\n');
fclose(fid);

disp('Done!');
tElapsed = toc(tStart);

fprintf('\nComputation completed in:\n');
fprintf('  %.2f seconds\n', tElapsed);
fprintf('  %.2f minutes\n', tElapsed/60);
fprintf('  %.2f hours\n', tElapsed/3600);

%% 11. Plot --------------------------------------------------------------
% Plotting 3D microstructure SURFACE
[Nx,Ny,Nz] = size(phase);

surface_mask = false(Nx,Ny,Nz);
surface_mask(1,:,:)   = true;
surface_mask(Nx,:,:)  = true;
surface_mask(:,1,:)   = true;
surface_mask(:,Ny,:)  = true;
surface_mask(:,:,1)   = true;
surface_mask(:,:,Nz)  = true;

phase_flat = phase(:);
X_flat = X(:); Y_flat = Y(:); Z_flat = Z(:);
surface_inds = find(surface_mask(:));
surface_inds = surface_inds(1:stride:end);
grain_inds = phase_flat(surface_inds);

valid = grain_inds >= 1 & grain_inds <= size(cmap_ipf,1);
surface_inds = surface_inds(valid);
grain_inds = grain_inds(valid);
pt_colors = cmap_ipf(grain_inds, :);

figure; clf; hold on;
scatter3(X_flat(surface_inds), Y_flat(surface_inds), Z_flat(surface_inds), 9, ...
    pt_colors, 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0);

if exist('n_pores','var') && n_pores > 0
    [xx,yy,zz]=sphere(30);
    for i=1:n_pores
        surf(x(i)+r(i)*xx, y(i)+r(i)*yy, z(i)+r(i)*zz, ...
            'FaceAlpha',0.7,'EdgeColor','none', ...
            'FaceColor',[0.4,0.4,0.6]);
    end
end

hx = xlabel('X_1 (mm)'); set(hx, 'Rotation', 21);
hy = ylabel('X_2 (mm)'); set(hy, 'Rotation', -21);
zlabel('X_3 (mm)');
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

%% Selective Isosurface Plot (Grains + True Pore Shapes)
pore_color     = [0.15, 0.15, 0.15];
pore_id        = N_grains + 1;

figure; clf; hold on;
n_show = min(num_grains_iso, N_grains);
cmap_grains = lines(max(64, n_show));
if n_show > size(cmap_grains,1)
    cmap_grains = parula(n_show);
end

for k = 1:n_show
    grain_mask = (phase == k);

    if any(grain_mask(:))
        fv = isosurface(xvec, yvec, zvec, ...
            smooth3(double(grain_mask),'box',3), 0.5);
        if ~isempty(fv.vertices)
            fv_red = reducepatch(fv, 0.15);
            patch(fv_red, 'FaceColor', cmap_grains(k,:), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.38, ...
                'SpecularStrength', 0.2);
        end
    end
end

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
        patch(fv_pores, 'FaceColor', pore_color, 'EdgeColor', 'none', ...
            'FaceAlpha', 1.0, 'DiffuseStrength', 0.8, ...
            'SpecularStrength', 0.5, 'SpecularExponent', 10);
    end
else
    n_plot_actual = 0;
    fprintf('No pores found in phase matrix.\n');
end

view([-42.71 22.97]);
axis equal tight; grid on; box on
camlight('headlight');
lighting gouraud;
material dull;
hx = xlabel('X_1 (mm)'); set(hx, 'Rotation', 21);
hy = ylabel('X_2 (mm)'); set(hy, 'Rotation', -21);
zlabel('X_3 (mm)');
xticks(linspace(xmin, xmax, 6));
yticks(linspace(ymin, ymax, 6));
zticks(linspace(zmin, zmax, 6));
xlim([xmin,xmax]);
ylim([ymin,ymax]);
zlim([zmin,zmax]);
set(gca, 'FontSize', 15, 'LineWidth', 1.3);
hAxes = findobj(gcf,"Type","axes");
hAxes.BoxStyle = "full";
hold off;


%% Pole figure
N = numel(ori_vec) / 3;
K = min(K_grain, N);
idx_grains = randperm(N, K);
ori_idx = bsxfun(@plus, 3*(idx_grains'-1), (1:3));
ori_idx = ori_idx';
ori_idx = ori_idx(:);
ori_vec_downsampled = ori_vec(ori_idx);
ori_vec_matrix = reshape(ori_vec_downsampled, 3, []).';
ori = orientation.byRodrigues(ori_vec_matrix, cs, ss);
psi = SO3vonMisesFisherKernel('halfwidth',5*degree);
odf = unimodalODF(ori, psi);
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

figure; clf;
plotPDF(odf, [Miller(1,0,0,cs), Miller(1,1,0,cs), Miller(1,1,1,cs)], ...
    'antipodal', 'contourf', 'smooth', 'resolution', 1*degree, 'minmax');
colorbar;

%% ---- END ----