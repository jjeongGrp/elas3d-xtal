# ELAS3D-Xtal

**An OpenMP-accelerated crystal elasticity solver with automated experiment-driven microstructure generation.**

ELAS3D-Xtal is a high-performance computational package for simulating the mechanical behavior of additively manufactured synthetic three-dimensional microstructures. It integrates a custom MATLAB-based microstructure generator with a modernized, parallelized version of the NIST finite element solver (originally NISTIR 6269).


---

## 📂 Repository Structure

```text
elas3d-xtal/
├── src/
│   ├── elas3dxtal_pcg.f90          # Main PCG Solver (OpenMP + HDF5)
│   └── microstructure_gen.m        # MATLAB Microstructure Generator
│
├── applications/
│   ├── sec3_1/                     # Section 3.1: Performance Scaling
│   │   ├── serial/                 # Serial Benchmarks
│   │   │   ├── n_100x100x100/
│   │   │   ├── n_200x200x200/
│   │   │   ├── n_300x300x300/
│   │   │   ├── n_400x400x400/
│   │   │   └── n_500x500x500/
│   │   ├── cg/                     # CG (OpenMP) Benchmarks
│   │   │   ├── n_100x100x100/
│   │   │   ├── n_200x200x200/
│   │   │   ├── n_300x300x300/
│   │   │   ├── n_400x400x400/
│   │   │   └── n_500x500x500/
│   │   └── pcg/                    # Preconditioned CG (OpenMP) Benchmarks
│   │       ├── n_100x100x100/
│   │       ├── n_200x200x200/
│   │       ├── n_300x300x300/
│   │       ├── n_400x400x400/
│   │       └── n_500x500x500/
│   │
│   ├── sec3_2/                     # Section 3.2: Anisotropy Validation
│   │   ├── inp/                    # Input templates
│   │   ├── A_0_5/                  # Anisotropy Ratio A = 0.5
│   │   ├── A_0_75/                 # Anisotropy Ratio A = 0.75
│   │   ├── A_1_0/                  # Anisotropy Ratio A = 1.0 (Isotropic)
│   │   ├── A_1_5/                  # Anisotropy Ratio A = 1.5
│   │   ├── A_2_0/                  # Anisotropy Ratio A = 2.0
│   │   ├── A_3_0/                  # Anisotropy Ratio A = 3.0
│   │   ├── A_4_0/                  # Anisotropy Ratio A = 4.0
│   │   └── A_5_0/                  # Anisotropy Ratio A = 5.0
│   │
│   └── sec3_3/                     # Section 3.3: AM Microstructure Simulation
│       ├── sim_EOS/                # Simulation: Gas Porosity
│       ├── sim_LoF/                # Simulation: Lack of Fusion Defects
│       ├── sim_KH/                 # Simulation: Keyhole Defects
│       ├── EOS_micro/              # Microstructure Data: EOS
│       ├── LoF_micro/              # Microstructure Data: Lack of Fusion
│       └── KH_micro/               # Microstructure Data: Keyhole
│
├── LICENSE
└── README.md

```

---

## 🛠️ Prerequisites (Software Requirements)

To compile and run ELAS3D-Xtal on Windows, please install the following software in the order listed below:

### Visual Studio 2022 (Community, Pro, or Enterprise)
* **Required For:** The C++ linker (`link.exe`) and build tools required by the Intel Fortran Compiler.
* **Download:** [Visual Studio 2022](https://visualstudio.microsoft.com/downloads/)
* **Installation Note:** During installation, select the following **Workloads**:
    *  **Desktop development with C++**
    * ️ **.NET desktop development**
    *️  **Python development**

### Intel Fortran Compiler (`ifx`)
* **Required For:** Compiling the Fortran source code with OpenMP support.
* **Download:** [Intel oneAPI HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)
* **Installation Note:** i. Run the installer and select **Custom Installation**.
    ii. When prompted *"Would you like to integrate with an IDE?"*, ensure the checkbox for **Microsoft Visual Studio 2022** is selected.

### HDF5 Library (Version 2.0.0)
* **Required For:** High-performance binary I/O (handling large voxel grids).
* **Download:** [HDF5 2.0.0 Release on GitHub](https://github.com/HDFGroup/hdf5/releases/tag/hdf5_2.0.0)
* **Specific File:** `hdf5-2.0.0-win-vs2022_intel.msi`
* **SHA256 Hash:** `8c625b68cb9b429208391f33b4ae675513c65c226ad3fe2a10226d93b16e2d35`
* **Installation Path:** The build commands assume the default location: `C:\Program Files\HDF_Group\HDF5\2.0.0\`

### MATLAB
* **Required For:** Generating microstructure inputs (`.h5` files) using the `microstructure_gen.m` script.
* **Recommended Version:** R2021b or newer.

**A. Required MathWorks Toolboxes:**
You must have the following toolboxes installed (check via `ver` command in MATLAB):
* **Statistics and Machine Learning Toolbox** (required for `lognrnd` and grain statistics).
* **Image Processing Toolbox** (required for voxel manipulation and segmentation).

**B. External Libraries (MTEX):**
The code relies on **MTEX** for crystallographic texture analysis and orientation mapping.
1.  **Download:** [MTEX Toolbox Website](https://mtex-toolbox.github.io/)
2.  **Install:**
    * Extract the folder (e.g., to `C:\Matlab_Toolboxes\mtex`).
    * Open MATLAB and run:
        ```matlab
        cd 'C:\Matlab_Toolboxes\mtex'
        startup_mtex
        ```
    * Run `microstructure_gen.m` only *after* MTEX is initialized.

---

## ⚙️ Compilation & Usage

### 1. Compilation (Windows)
Open your **Intel oneAPI Command Prompt for Intel 64** (search for it in the Start Menu—do **not** use the standard `cmd.exe`) and navigate to the `src` directory. Run **one** of the following commands to build the solver.

**A. PCG Solver (Recommended)**
High-performance OpenMP solver with point-block Jacobi preconditioning.

```cmd
ifx /O3 /QxHost /Qunroll /Qopenmp /Qipo /I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" elas3dxtal_pcg.f90 /Fe:elas3dxtal_pcg.exe /link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" hdf5_fortran.lib hdf5.lib

```

**B. CG Solver**
Standard Conjugate Gradient solver without preconditioning.

```cmd
ifx /O3 /QxHost /Qunroll /Qopenmp /Qipo /I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" elas3dxtal_cg.f90 /Fe:elas3dxtal_cg.exe /link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" hdf5_fortran.lib hdf5.lib

```

**C. Serial Solver (Validation)**
Single-core baseline for debugging or small grids.

```cmd
ifx /I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" elas3dxtal_serial.f90 /Fe:elas3dxtal_serial.exe /link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" hdf5_fortran.lib hdf5.lib

```


### 2. Execution Guide

**Step A: Generate Microstructure (MATLAB)**

1. Open MATLAB and navigate to the `src/` folder.
2. Ensure the **MTEX** toolbox is active (`startup_mtex`).
3. Run the generation script:
```matlab
elas3dxtal_input

```

* **Action:** Generates a 3D voxel mesh based on grain statistics and pore data.
* **Output:** Creates `input_structure_poly.h5` (required by the solver).



**Step B: Run the Solver (Command Line)**
In your Intel Command Prompt, execute the compiled solver:

```cmd
elas3dxtal_pcg.exe

```

* **Action:** Reads `input_structure_poly.h5` and solves the elastic equilibrium equations using OpenMP threads.
* **Output:** Generates `fullfield_poly.h5`, containing stress, phase, and orientation fields.


**Step C: Post-Processing (MATLAB)**

1. In MATLAB, run the post-processing script:
```matlab
elas3dxtal_postprocessing

```


* **Action:** Reads `fullfield_poly.h5`.
* **Output:** Visualizes 3D von Mises stress maps and cross-sectional contour plots.



---

## 🔧 Configuration Details: Linking MATLAB to Fortran

To run a simulation, ensure that the parameters in the MATLAB generator (`microstructure_gen.m`) match the hardcoded parameters in the Fortran solver (`elas3dxtal_pcg.f90`).

### 1. `microstructure_gen.m` Configuration
This script generates the synthetic microstructure. The parameters in the "Configuration & Parameters" section control the output.

#### A. Domain & Grid Resolution
These variables define the physical size of the simulation volume and the voxel discretization.

* **`cube_len`**: The physical length of the cubic domain side (e.g., in mm).
* **`grid_num`**: The number of voxels along one axis ($N_x = N_y = N_z$).
    * *Note:* This determines the total problem size ($N^3$). Increasing this improves resolution but drastically increases memory usage and solve time.
* **`dx`**: The voxel size (automatically calculated as `cube_len / grid_num`).

#### B. Grain Statistics (Morphology)
Controls the size and shape of the grains using statistical distributions.

* **`mean_d`**: Average grain diameter (in physical units, e.g., mm).
    * *Note:* Smaller values result in more grains for a given `cube_len`.
* **`std_d`**: Standard deviation of the grain diameter (controls size variability).
* **`distribution_type`**: Statistical distribution used for grain sizes (typically `'lognormal'`).
* **`mean_aspect`**: Average grain aspect ratio (Length / Width).
    * *Value = 1.0:* Equiaxed grains.
    * *Value > 1.0:* Columnar/Elongated grains.
* **`std_aspect`**: Standard deviation of the aspect ratio.

#### C. Crystallographic Orientation
Defines how crystal orientations (texture) are assigned to the generated grains.

* **`orientation_type`**:
    * `'random'`: Assigns uniformly random orientations to all grains.
    * `'textured'`: Biases orientations towards a specific direction (mimics fiber texture).
* **`sigma_spread`**: Spread (in degrees) of the orientation distribution if `'textured'` is selected.

#### D. Defect (Pore) Insertion
Controls the inclusion of porosity or defects from experimental data.

* **`pore`**: Toggle (`'on'` / `'off'`) to include defects.
* **`input_pore_file`**: Path to the external file (Excel/CSV) containing defect data (Centers and Radii).
* **`pore_type`**: Descriptive label for the defect type (e.g., `'LoF'`, `'Keyhole'`).

#### E. Output Controls
Manages which files are generated for the solver and visualization.

* **`save_data`**: Toggle (`'on'`) to save the MATLAB workspace (`.mat` file).
* **`elas3d_input`**: **CRITICAL.** Must be set to `'on'` to generate the HDF5 file (`input_structure_poly.h5`) required by the Fortran solver.
* **`micro_plot`**: Toggle (`'on'`) to generate visualization figures (cross-sections, 3D renders) after generation.


### 2. Linking to `elas3dxtal_pcg.f90` (Fortran Solver)
After generating the microstructure in MATLAB and exporting the results to 'output_summary.txt', manually update the Fortran solver parameters to match your new geometry and material.

**CRITICAL CONCEPT: The "+1 Rule" for Grid Points**
* **MATLAB** defines the grid by **voxels** (elements).
* **Fortran** defines the grid by **nodes** (grid Points).
* *Relationship:* `Nodes = Voxels + 1`

#### A. Grid Resolution (Domain Size)
Set the array dimensions in the Fortran module to match the MATLAB grid.

* **`md`**: Total number of **nodes** along one axis.
    * *Formula:* `md = grid_num + 1`
    * *Example:* If MATLAB `grid_num = 100`, then Fortran `md = 101`.
* **`nx`, `ny`, `nz`**: Set all equal to `md` (for a cubic domain).

#### B. Phase & Grain Counts
The solver needs to know exactly how many distinct crystal grains exist to allocate memory for their orientations.

* **`n_grains`**: Total number of grains generated.
    * *Action:* Check the variable `n_grains` (or `num_grains`) in your MATLAB workspace after generation and copy that value here.
* **`nphase`**: The total number of phases in the system.
    * *Formula:* `nphase = n_grains + 1` (The `+1` accounts for the pore/void phase).

#### C. Material Stiffness (Constitutive Law)
MATLAB handles geometry; Fortran handles physics. Define the elastic properties of your specific material here.

* **`C11_local`, `C12_local`, `C44_local`**: The independent elastic stiffness constants for a cubic crystal.
    * *Units:* **Pascals (Pa)**.
    * *Note:* Ensure you convert from GPa to Pa (e.g., $200 \text{ GPa} = 200.0\text{d}9$).

#### D. Boundary Conditions (Applied Loading)
Controls the macroscopic strain applied to the simulation volume.

* **`aml`**: The magnitude of the applied load (e.g., `1.0e-3` for 0.1% strain).
* **`aml_exx`, `aml_eyy`, `aml_ezz`**: The specific components of the macroscopic strain tensor to enforce.
    * *Example:* To simulate uniaxial tension in Z, set `aml_ezz = aml` and others to `0.0d0`.



---

## 📊 Reproducing Paper Results

The `applications/` folder contains specific configuration files and scripts to reproduce the results presented in **Jeong & Sundararaghavan (2026)**.

	
### 1. Performance Scaling (Section 3.1)
* **Goal:** Validate ELAS3D-Xtal solver accuracy against the analytical Eshelby solution and benchmark the performance of different solvers (Serial, CG, PCG).
* **Location:** `applications/sec3_1/`
* **Experiments:**
    i.  **Accuracy Validation:**
        * **Setup:** Single spherical inclusion in an isotropic matrix.
        * **Action:** Compare the numerical von Mises stress field ($\sigma_{vm}$) along the inclusion centerline against the exact **analytical Eshelby Green's function solution**.
    ii.  **Performance Scaling:**
        * **Solvers:** Compare convergence and execution time for:
            * `elas3dxtal_serial.exe` (Baseline)
            * `elas3dxtal_cg.exe` (Standard Conjugate Gradient)
            * `elas3dxtal_pcg.exe` (Preconditioned CG)
        * **Grid Sizes:** Run simulations on $100^3$, $200^3$, $300^3$, $400^3$, and $500^3$ voxel grids to quantify $O(N \log N)$ scaling.
		
### 2. Anisotropy Validation (Section 3.2)
* **Goal:** Verify that the ELAS3D-Xtal solver correctly handles high degrees of elastic anisotropy, and investigate how the Zener anisotropy ratio ($A$) alters stress concentrations ($K_{t}$) near defects in a polycrystalline matrix with random texture. 
* **Location:** `applications/sec3_2/`
* **Configuration:**
    * **Variable:** Zener Anisotropy Ratio, $A = \frac{2C_{44}}{C_{11} - C_{12}}$.
    * **Test Cases:** Subdirectories for $A = \{0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0\}$.
* **Action:**
    i.  **Compile:** Build the solver using the specific stiffness constants ($C_{11}, C_{12}, C_{44}$) provided in the `elas3dxtal_pcg.f90` of each subdirectory.
    ii.  **Run:** Execute `elas3dxtal_pcg.exe` on the single-inclusion benchmark.
    iii.  **Compare:** Plot the resulting von Mises stress field $\sigma_{vm}$ to observe the directionality of stress concentrations in highly anisotropic medium.
	
	
### 3. Application to AM Microstructures (Section 3.3)
* **Goal:** Demonstrate full-field stress simulation on realistic additively manufactured (AM) microstructures containing process-induced defects.
* **Location:** `applications/sec3_3/`
* **Test Cases (Microstructure Generation):**
    1.  **EOS (Gas Porosity):** `eos_micro/` - Microstructure containing small, spherical voids and exhibiting strong **[101] fiber texture**.
    2.  **LoF (Lack of Fusion):** `lof_micro/` - Microstructure containing large, flat, irregular (oblate) voids  and **random texture** (no preferred orientation).
    3.  **Keyhole:** `keyhole_micro/` - Microstructure containing deep, narrow (prolate) voids  and exhibiting strong **[101] fiber texture**.
* **Workflow:**
    1.  **Generate:** Run `elas3dxtal_input.m` in the respective `sim_*/` folder to generate the voxel mesh from experimental defect data.
    2.  **Simulate:** Run `elas3dxtal_pcg.exe` to solve for elastic equilibrium.
    3.  **Visualize:** Run `elas3dxtal_postprocessing.m` to compute the von Mises stress field $\sigma_{vm}$ and visualize the "hotspots" near pore boundaries.

---

## 📄 Citation

If you use this code in your research, please cite **both** the new methodology and the original NIST algorithm:

**1. ELAS3D-Xtal (This Work):**

> Jeong, J., & Sundararaghavan, V. (2026). *ELAS3D-Xtal: An OpenMP-accelerated crystal elasticity solver with automated experiment-driven microstructure generation*. arXiv preprint arXiv:2602.07354. https://doi.org/10.48550/arXiv.2602.07354

**2. Original NIST Solver:**

> Garboczi, E. (1998). *Finite Element and Finite Difference Programs for Computing the Linear Electric and Elastic Properties of Digital Images of Random Materials*. NIST Interagency/Internal Report (NISTIR) 6269.

### BibTeX

```bibtex
@article{jeong2026elas3d,
  title={ELAS3D-Xtal: An OpenMP-accelerated crystal elasticity solver with automated experiment-driven microstructure generation},
  author={Jeong, Juyoung and Sundararaghavan, Veera},
  journal={arXiv preprint arXiv:2602.07354},
  year={2026},
  doi={10.48550/arXiv.2602.07354},
  url={https://doi.org/10.48550/arXiv.2602.07354}
}

@article{garboczi1998finite,
  title={Finite element and finite difference programs for computing the linear electric and elastic properties of digital images of random materials},
  author={Garboczi, Edward J},
  year={1998},
  journal={Report NISTIR},
  publisher={Edward J. Garboczi}
}

```

---

## ⚖️ License & Attribution

This project uses a **hybrid license** model to respect the original government work:

1. **New contributions:** Licensed under the **MIT License**.
2. **Original Solver:** Derived from NISTIR 6269. This software is a derivative work of NISTIR 6269, republished courtesy of the National Institute of Standards and Technology.

**Full details:** See the [`LICENSE`] file.

---

## 🤝 Acknowledgements

This work was supported by the Defense Advanced Research Projects Agency (DARPA) SURGE program under Cooperative Agreement No. HR0011-25-2-0009, "Predictive Real-time Intelligence for Metallic Endurance (PRIME)."

*Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of DARPA or the United States Government.*

---

## 📫 Contact

For questions regarding the code, methodology, or the associated publication, please feel free to reach out:

* **Juyoung Jeong** - jjuyoung@umich.edu
* **Prof. Veera Sundararaghavan** - veeras@umich.edu
* **Affiliation:** Department of Aerospace Engineering, University of Michigan

**Bug Reports & Feature Requests:** If you encounter any issues while compiling or running the solver, please use the [GitHub Issues](https://github.com/jjeongGrp/elas3d-xtal/issues) page to report them.

**Project Link:** [https://github.com/jjeongGrp/elas3d-xtal](https://github.com/jjeongGrp/elas3d-xtal)
