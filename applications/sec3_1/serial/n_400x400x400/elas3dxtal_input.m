%% Eshelby Microstructure Generator for ELAS3D-Xtal
% ------------------------------------------------------------------------------
% This script generates a digital 3D voxelized domain containing a single 
% ellipsoidal inclusion embedded within a matrix. It serves as a pre-processing 
% module for computational solid mechanics solvers, ELAS3D-Xtal, to 
% verify numerical results against exact analytical solutions of the classic 
% Eshelby inclusion problem.
%
% Dependency:
%   - Base MATLAB (No external toolboxes required)
%
% Inputs:
%   - Physical domain size and grid resolution parameters.
%   - Inclusion geometry parameters (semi-axes and aspect ratio).
%
% Outputs:
%   - 'input_structure_poly.h5': HDF5 file containing the 3D phase map and 
%     default crystallographic orientations for solver ingestion.
%   - 'params_Eshelby_ellipsoidal.mat': Workspace data containing geometric bounds.
%   - Visualization figures: 2D mid-plane cross-section, 3D smooth isosurface 
%     render, and 3D voxelized scatter plot.
%
% Core Workflow:
%   1. Define computational domain and voxel discretization.
%   2. Define inclusion parameters and generate full 3D coordinate grids.
%   3. Apply an analytical mathematical mask for the ellipsoid to assign 
%      phase IDs (Matrix = 1, Inclusion = 2).
%   4. Export flattened phase array and default orientations to HDF5 format.
%   5. Generate 2D and 3D visual diagnostics of the generated microstructure.
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

% =========================================================================
%% 1. Physical Domain Size and Mesh Resolution
% =========================================================================
% IMPORTANT NOTE ON "n", VOXELS, AND GRID POINTS (MATCHING THE ANALYTICAL CONVENTION)
% ------------------------------------------------------------------------------
% For strict voxel meshes in 1D:
%   - n voxels (cells) span the domain with n intervals
%   - which requires n+1 edge nodes
%
% HOWEVER, the analytical Eshelby *reference solution used for verification*
% in this workflow samples fields on an "n-point linspace grid" that INCLUDES
% both endpoints of the domain (i.e., it treats the domain as n *sample points*,
% not n voxels with n+1 nodes).
%
% To ensure the numerical microstructure is generated on the same physical
% sampling locations (and avoids a half-cell shift when comparing to the
% analytical reference), this generator follows the SAME sampling convention:
%   - xvec/yvec/zvec contain EXACTLY n points spanning [-L, +L] (endpoints included)
%   - implied spacing is dx = 2L/(n-1)
%
% If the solver expects the strict "n voxels -> n+1 nodes" convention, 
% construct edge grids (n+1) and evaluate the phase at voxel centers (n). 
% That is intentionally NOT used here, to remain consistent with
% the analytical reference sampling convention.
% ------------------------------------------------------------------------------
L = 50.0;         % Half-width of domain in all directions (um)
n = 100;          % Number of GRID SAMPLE POINTS per direction (endpoints included).
                  % NOTE: This is *not* "n voxels with n+1 nodes".
                  % Implied interval count is (n-1), and dx = 2L/(n-1).

nx = n;           % Number of sample points in x
ny = n;           % Number of sample points in y
nz = n;           % Number of sample points in z

% =========================================================================
%% 2. Inclusion Parameters (Ellipsoid Geometry)
% =========================================================================
alp = 0.999;      % Aspect ratio a1/a2 (1 = sphere, <1 = oblate, >1 = prolate)
a2 = 2;           % Second semi-axis (arbitrary units, must match L scale)
a1 = a2 * alp;    % First semi-axis (calculated via aspect ratio)
a3 = a2;          % Third semi-axis (assumed symmetric with a2 here)
nphase = 2;       % Total phases: 1 = matrix, 2 = inclusion

% Generate 1D coordinate vectors (n sample points INCLUDING endpoints)
% This matches the analytical reference-solution sampling convention.
xvec = linspace(-L, L, nx);   
yvec = linspace(-L, L, ny);   
zvec = linspace(-L, L, nz);   

% Generate full 3D coordinate grids at these sample points
[X1, X2, X3] = ndgrid(xvec, yvec, zvec);

% Inclusion mask using the standard ellipsoid equation:
% (x1/a1)^2 + (x2/a2)^2 + (x3/a3)^2 <= 1  ==> inside inclusion
incl_mask3d = (X1/a1).^2 + (X2/a2).^2 + (X3/a3).^2 <= 1;

% Initialize phase array (Default: Matrix = phase 1)
phase = ones(nx, ny, nz, 'int32');

% Assign the inclusion ID to the masked region (Inclusion = phase 2)
phase(incl_mask3d) = 2;



% =========================================================================
%% 3. Export Phase and Orientation to HDF5
% =========================================================================
h5filename = 'input_structure_poly.h5';

% Overwrite protection
if isfile(h5filename)
    delete(h5filename);
    fprintf('Previous file "%s" deleted.\n', h5filename);
end

% Flatten 3D phase array to 1D column vector for HDF5
pix_flat = int32(phase(:));
h5create(h5filename, '/pix', [length(pix_flat) 1], 'Datatype', 'int32');
h5write(h5filename, '/pix', pix_flat);
h5writeatt(h5filename, '/pix', 'dimensions', int32([nx ny nz]));

% Handle crystallographic orientations
% If an orientation vector is provided externally, write it.
% Otherwise, write zeros (default orientation) for the non-matrix phase(s).
if exist('ori_vec', 'var')
    ovec = ori_vec;
else
    ovec = zeros(3*(nphase-1), 1);
end

h5create(h5filename, '/orientation', size(ovec), 'Datatype', 'double');
h5write(h5filename, '/orientation', ovec);
h5writeatt(h5filename, '/orientation', 'dimensions', int32(size(ovec)));

disp(['Saved HDF5 input structure to ', h5filename]);

% Save basic geometry parameters to MATLAB workspace file
save('params_Eshelby_ellipsoidal.mat', 'nx', 'ny', 'nz', 'L', ...
    'xvec', 'yvec', 'zvec', 'a1', 'a2', 'a3');

% =========================================================================
%% 4. Microstructure Plotting
% =========================================================================

% --- FIGURE 1: 2D Cross-Section (Mid-plane slice) ---
figure; clf; 
center_x = round(nx/2);
% Extract slice at X_mid, transpose for standard image viewing
ph2d = squeeze(phase(center_x, :, :))';

imagesc(yvec, zvec, ph2d);  
colormap([[1 1 1]; [0 0 1]]);  % Matrix = White, Inclusion = Blue 
axis image;
clim([1 2]);

axis_range = [-L, L];
interval = L/2;
xlim(axis_range); ylim(axis_range); 
xticks(axis_range(1):interval:axis_range(2));
yticks(axis_range(1):interval:axis_range(2));

xlabel('Position X_2 (μm)');
ylabel('Position X_3 (μm)');
set(gca,'FontSize',15)

% Draw a faint grid to represent the finite element mesh
hold on
nyy = length(yvec); nzz = length(zvec);
for k = 1:nyy
    plot([yvec(1) yvec(end)], zvec([k k]), 'k-', 'LineWidth', 0.01);
end
for j = 1:nzz
    plot(yvec([j j]), [zvec(1) zvec(end)], 'k-', 'LineWidth', 0.01);
end
hold off

% --- FIGURE 2: 3D Smooth Isosurface Render ---
figure; clf;
[X, Y, Z] = ndgrid(xvec, yvec, zvec);
inclu = (phase == 2);

% Render the inclusion boundary smoothly
p = patch(isosurface(X, Y, Z, inclu, 0.5));
set(p, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.8);

xlabel('X_1-axis (\mum)', 'Rotation',30, 'FontWeight', 'bold');
ylabel('X_2-axis (\mum)', 'Rotation',-30, 'FontWeight', 'bold');
zlabel('X_3-axis (\mum)', 'FontWeight', 'bold');
set(gca,"FontSize",15)
axis equal tight
view([-41.09 35.08]);
camlight; lighting gouraud
hold on

% Render a translucent bounding box for the matrix domain
Vtx = [xvec(1) yvec(1) zvec(1);
       xvec(end) yvec(1) zvec(1);
       xvec(end) yvec(end) zvec(1);
       xvec(1) yvec(end) zvec(1);
       xvec(1) yvec(1) zvec(end);
       xvec(end) yvec(1) zvec(end);
       xvec(end) yvec(end) zvec(end);
       xvec(1) yvec(end) zvec(end)];
CubeFaces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
patch('Vertices',Vtx,'Faces',CubeFaces, 'FaceColor','w','FaceAlpha',0.15, ...
    'EdgeColor','none');

plot3([xvec(1) xvec(end) xvec(end) xvec(1) xvec(1)], ...
      [yvec(1) yvec(1) yvec(end) yvec(end) yvec(1)], ...
      [zvec(1) zvec(1) zvec(1)  zvec(1)  zvec(1)], 'k-', 'LineWidth', 1.0);
plot3([xvec(1) xvec(end) xvec(end) xvec(1) xvec(1)], ...
      [yvec(1) yvec(1) yvec(end) yvec(end) yvec(1)], ...
      [zvec(end) zvec(end) zvec(end) zvec(end) zvec(end)], ...
      'k-', 'LineWidth', 1.0); % top
for xnow = [xvec(1), xvec(end)]
    for ynow = [yvec(1), yvec(end)]
        plot3([xnow xnow],[ynow ynow],[zvec(1) zvec(end)], ...
            'k-', 'LineWidth', 1.0)
    end
end
hold off
box on; grid on;
legend('Inclusion','Matrix','Location','northeast')


% --- FIGURE 3: 3D Voxel Render (Scatter) ---
figure; clf;
inclu = (phase == 2);
[xg, yg, zg] = ndgrid(xvec, yvec, zvec);
idx = find(inclu);

scatter3(xg(idx), yg(idx), zg(idx), 15, [0 0 1], 's', 'filled');
hold on

plot3([xvec(1) xvec(end) xvec(end) xvec(1) xvec(1)], ...
      [yvec(1) yvec(1) yvec(end) yvec(end) yvec(1)], ...
      [zvec(1) zvec(1) zvec(1)  zvec(1)  zvec(1)], 'k-', 'LineWidth', 1.5);
plot3([xvec(1) xvec(end) xvec(end) xvec(1) xvec(1)], ...
      [yvec(1) yvec(1) yvec(end) yvec(end) yvec(1)], ...
      [zvec(end) zvec(end) zvec(end) zvec(end) zvec(end)], ...
      'k-', 'LineWidth', 1.5);
for xnow = [xvec(1), xvec(end)]
    for ynow = [yvec(1), yvec(end)]
        plot3([xnow xnow],[ynow ynow],[zvec(1) zvec(end)],'k-', ...
            'LineWidth', 1.5)
    end
end
n = length(xvec);   
step = max(floor(n/20),1);

for k = 1:step:n
    plot3([xvec(1) xvec(1)], [yvec(k) yvec(k)], [zvec(1) zvec(end)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(1) xvec(1)], [yvec(1) yvec(end)], [zvec(k) zvec(k)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(end) xvec(end)], [yvec(k) yvec(k)], [zvec(1) zvec(end)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(end) xvec(end)], [yvec(1) yvec(end)], [zvec(k) zvec(k)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);

    plot3([xvec(k) xvec(k)], [yvec(1) yvec(1)], [zvec(1) zvec(end)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(1) xvec(end)], [yvec(1) yvec(1)], [zvec(k) zvec(k)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(k) xvec(k)], [yvec(end) yvec(end)], [zvec(1) zvec(end)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(1) xvec(end)], [yvec(end) yvec(end)], [zvec(k) zvec(k)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);

    plot3([xvec(k) xvec(k)], [yvec(1) yvec(end)], [zvec(1) zvec(1)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(1) xvec(end)], [yvec(k) yvec(k)], [zvec(1) zvec(1)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(k) xvec(k)], [yvec(1) yvec(end)], [zvec(end) zvec(end)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
    plot3([xvec(1) xvec(end)], [yvec(k) yvec(k)], [zvec(end) zvec(end)], ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.05);
end

hold off
xlabel('X_1-axis (μm)', 'Rotation',30, 'FontWeight', 'bold');
ylabel('X_2-axis (μm)', 'Rotation',-30, 'FontWeight', 'bold');
zlabel('X_3-axis (μm)', 'FontWeight', 'bold');
set(gca,"FontSize",15)
axis equal tight
view([-41.09 35.08]);
box on; grid on;
