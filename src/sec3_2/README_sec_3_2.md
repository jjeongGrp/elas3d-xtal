# Section 3.2: Anisotropy Validation

## Execution Workflow

### Step 1: Generate the Base Microstructure

Generate the single master microstructure in the `inp/` directory.

1. Open MATLAB and navigate to `applications/sec3_2/inp/`.
2. Run the microstructure generator:
```matlab
elas3dxtal_input

```

3. Verify that `input_structure_poly.h5` has been successfully generated in the `inp/` folder.

### Step 2: Distribute the Input File

Copy the newly generated `input_structure_poly.h5` and `xct_poly_params_SS316L.mat` from the `inp/` folder and paste it into **every** target simulation folder:

* `A_0_5/`
* `A_0_75/`
* `A_1_0/`
* `A_1_5/`
* `A_2_0/`
* `A_3_0/`
* `A_4_0/`
* `A_5_0/`


### Step 3: Compile and Run the Solvers

For each `A_*` folder, compile its local Fortran source file, as the stiffness constants are hardcoded into each specific file.

1. Open the **Intel oneAPI Command Prompt**.
2. Navigate into a specific folder (e.g., `cd A_0_5`).
3. **Compile** the code:
```cmd
ifx /O3 /QxHost /Qunroll /Qopenmp /Qipo /I"C:\Program Files\HDF_Group\HDF5\2.0.0\mod\shared" elas3dxtal_pcg.f90 /Fe:elas3dxtal_pcg.exe /link /LIBPATH:"C:\Program Files\HDF_Group\HDF5\2.0.0\lib" hdf5_fortran.lib hdf5.lib

```


4. **Run** the executable:
```cmd
elas3dxtal_pcg.exe

```


5. Repeat this process for all `A_*` folders.

### Step 4: Post-Processing & Visualization

Once the solver finishes, it outputs `fullfield_poly.h5`.

1. Open MATLAB and navigate into the specific folder (e.g., `A_0_5/`).
2. Run the post-processing script:
```matlab
elas3dxtal_postprocessing

```



```