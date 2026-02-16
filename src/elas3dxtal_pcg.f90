!--------------------------------------------------------------------------------
! ELAS3D-Xtal: Preconditioned Conjugate Gradient Solver for 3D Crystal Elasticity
! -------------------------------------------------------------------------------
! This program solves the 3D linear elastic response of a voxelized multi-phase
! microstructure (e.g., XCT-derived defect-containing metal, polycrystals) using 
! a regular grid finite element method. The code supports both isotropic and 
! anisotropic phases and applies periodic boundary conditions. The stress and
! strain response to user-prescribed macroscopic strains is computed via a 
! preconditioned conjugate gradient (PCG) energy minimization algorithm, heavily
! optimized for shared-memory parallelization using OpenMP.
!    
! This program builds upon the foundational algorithms of:
!   Garboczi, E. J. (1998). Finite element and finite difference programs for 
!   computing the linear electric and elastic properties of digital images of 
!   random materials. NISTIR 6269.    
!
! Major Capabilities:
!   - Reads large-scale microstructure (phase/texture) data via HDF5.
!   - Assembles periodic 3D finite element stiffness matrices.
!   - Supports isotropic matrix/inclusion and highly anisotropic polycrystals.
!   - Applies macroscopic loading (strain-controlled).
!   - Solves for equilibrium using an OpenMP-accelerated PCG method.
!   - Extracts and writes full-field stress and strain tensors for all voxels.
!   - Outputs solver history, effective moduli, and phase statistics.
!
! Authors:       Juyoung Jeong and Veera Sundararaghavan
! Affiliation:   Department of Aerospace Engineering, University of Michigan,
!                Ann Arbor, MI 48109, USA.
! Repository:    https://github.com/jjeongGrp/elas3d-xtal
!
! Acknowledgement:
!   This work was supported by the Defense Advanced Research Projects Agency 
!   (DARPA) SURGE program under Cooperative Agreement No. HR0011-25-2-0009, 
!   'Predictive Real-time Intelligence for Metallic Endurance (PRIME)'.
!
! License:       MIT License (See LICENSE file in the repository)
!
! Revision History:
!   July 28, 2025      - Initial version.
!   August 22, 2025    - Memory optimizations and HDF5 I/O integration.
!   August 28, 2025    - Chunking iteration (dembx_OpenMP).
!   November 15, 2025  - Memory access pattern optimizations.
!                      - Load balancing via OpenMP guided scheduling.
!                      - Implementation of the Preconditioned Conjugate 
!                        Gradient (PCG) algorithm.
!------------------------------------------------------------------------------  
      module elas3d_mod
      ! ==============================================================================
      ! MODULE: elas3d_mod
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Defines the core data structures and shared variables for the 3D linear 
      !   elasticity finite element solver. By encapsulating grid parameters, 
      !   material properties, and the large arrays required for the FEM assembly 
      !   into a custom type (elas3d_data_type), this module ensures thread-safe 
      !   data passing of input and output variables, eliminating the need for 
      !   implicit global variables.
      !
      ! DATA STRUCTURE COMPONENTS (elas3d_data_type):
      !   - pix (intent(in)): 1D array representing the 3D voxel phase map.
      !   - ib (intent(in)): Periodic neighbor connectivity mapping array.
      !   - exx, eyy, ezz, exy, exz, eyz (intent(in)): Applied macroscopic strains.
      !   - rotatedstiffness (intent(in)): 1D array of pre-calculated rotated grain tensors.
      !   - cmod (intent(out)): 6x6 local constitutive stiffness matrices for each phase.
      !   - dk (intent(out)): Integrated local finite element stiffness matrices.
      !   - b (intent(out)): Global right-hand-side force vector for the solver.
      !   - C (intent(out)): Global energy offset scalar representing base strain energy.
      ! ==============================================================================	  
       use, intrinsic :: iso_fortran_env, dp=>real64 
       implicit none
       
       ! ======================================================================
       ! ---- 1. Domain & Grid Resolution -------------------------------------
       ! Defines the finite element grid. 
       ! Note: Fortran defines nodes, MATLAB defines voxels. Nodes = Voxels + 1.
       ! ======================================================================
       integer, parameter :: md = 101     ! Number of nodes per direction (e.g., 100 voxels + 1)
       integer, parameter :: nx = md      ! Number of grid points in X-direction
       integer, parameter :: ny = md      ! Number of grid points in Y-direction
       integer, parameter :: nz = md      ! Number of grid points in Z-direction
       integer, parameter :: ns = nx * ny * nz   ! Total number of nodal sites in the domain

       ! ======================================================================
       ! ---- 2. Microstructure Definition ------------------------------------
       ! Controls phase allocation and FE element properties.
       ! ======================================================================
       integer, parameter :: n_grains = 910      ! Total number of distinct crystal grains
       integer, parameter :: nphase = n_grains+1 ! Total phases (Grains + 1 Inclusion/Pore phase)
       integer, parameter :: nphmax = nphase     ! Maximum number of phases (used for array bounds)   
       integer, parameter :: nfaces = 27         ! Number of periodic neighbor interactions (3x3x3 grid) 
       integer, parameter :: ndof = 3            ! Degrees of freedom per node (Displacements in X, Y, Z)
	   integer, parameter :: nnode_fe = 8        ! Nodes per finite element (8-node linear hexahedral/brick)
       integer, parameter :: ngauss   = 3        ! Gauss integration points per direction (3x3x3 = 27 total)       
       
       ! ======================================================================
       ! ---- 3. Material Properties (Single Crystal Stiffness) ---------------
       ! Independent elastic constants for a cubic crystal lattice (e.g., FCC).
       ! Must be defined in Pascals (Pa).
       ! ======================================================================
       real(dp), parameter :: C11_local = 206.00d9   ! C11 Stiffness component (Pa)
       real(dp), parameter :: C12_local = 133.00d9   ! C12 Stiffness component (Pa)
       real(dp), parameter :: C44_local = 119.00d9   ! C44 Stiffness component (Pa) 	   
	   
       ! ======================================================================
       ! ---- 4. Isotropic Equivalents & Phase Flags --------------------------
       ! Calculates Voigt-Reuss-Hill averages for preconditioner/baseline.
       ! ======================================================================
	   real(dp), parameter :: K0 = (C11_local + 2.0d0*C12_local)/3.0d0             ! Bulk Modulus 
       real(dp), parameter :: G0 = (C11_local - C12_local + 3.0d0*C44_local)/5.0d0 ! Shear Modulus
       real(dp), parameter :: E0 = 9.0d0*K0*G0/(3.0d0*K0 + G0)                     ! Young's Modulus, Matrix (Pa)  
       real(dp), parameter :: nu0 = (3.0d0*K0 - 2.0d0*G0) / (2.0d0*(3.0d0*K0+G0))  ! Poisson's Ratio, Matrix
	   
       ! Properties of the inclusion/defect phase (0.0 represents a void/pore)	   
       real(dp), parameter :: E1 = 0.0d0        ! Young's Modulus, Inclusion
       real(dp), parameter :: nu1 = 0.0d0       ! Poisson's Ratio, Inclusion      

       ! Solver flag for microstructure behavior	   
       integer, parameter :: flag_m = 1         ! 0: Treat matrix as single isotropic block
                                                ! 1: Treat matrix as textured polycrystal (rotates stiffness per grain)
       
       ! ======================================================================
       ! ---- 5. Boundary Conditions (Macroscopic Loading) --------------------
       ! Defines the far-field macroscopic strain applied to the periodic box.
       ! ======================================================================
       real(dp), parameter :: aml     =  1.0e-3  ! Applied load magnitude (e.g., 1.0e-3 = 0.1% strain)
       real(dp), parameter :: aml_exx =  0.00d0  ! Applied normal strain XX
       real(dp), parameter :: aml_eyy =  0.00d0  ! Applied normal strain YY    
       real(dp), parameter :: aml_ezz =  aml     ! Applied normal strain ZZ (Uniaxial tension in Z)
       real(dp), parameter :: aml_exz =  0.00d0  ! Applied shear strain XZ
       real(dp), parameter :: aml_eyz =  0.00d0  ! Applied shear strain YZ
       real(dp), parameter :: aml_exy =  0.00d0  ! Applied shear strain XY         
       
       ! ======================================================================
       ! ---- 6. PCG Solver Iteration Controls --------------------------------
       ! Configures the Preconditioned Conjugate Gradient minimization steps.
       ! ======================================================================
       ! Keep kmax as 1. For this specific linear elastic formulation, a single 
       ! continuous solve is optimal. Forcing global restarts clears the conjugate 
       ! search directions and severely slows down convergence without benefit.	   
       integer, parameter :: kmax = 1              ! Outer loop counter for global solver restarts
       integer, parameter :: ldemb = 10000         ! Max number of Conjugate Gradient iterations
       integer, parameter :: n_iter = kmax * ldemb ! Total allowable solver steps
       integer, parameter :: block_size = 8192     ! Memory chunk size for OpenMP cache optimization
       
       ! ======================================================================
       ! ---- 7. Convergence Criterion ----------------------------------------
       ! ======================================================================
       real, parameter :: tol = 1.0d-8             ! Energy gradient tolerance to determine equilibrium

       ! ======================================================================
       ! ---- 8. Global Data Structure (elas3d_data_type) ---------------------
       ! Groups all dynamically allocated arrays to pass data cleanly through
       ! subroutines without excessive argument lists.
       ! ======================================================================
       type :: elas3d_data_type
         ! FE/Solver Variables: Displacements (u), Forces (b), Gradients (gb, h)
         real(dp), allocatable  :: u(:,:), gb(:,:), b(:,:)
		 real(dp), allocatable  :: h(:,:)

         ! Stiffness arrays: local element moduli (cmod), global assembled moduli (dk)		 
         real(dp), allocatable  :: cmod(:,:,:), dk(:,:,:,:,:)
		 
         ! Nodal indexing and phase mappings		 
         integer, allocatable   :: ib(:,:), pix(:)
		 
         ! Macroscopic strain and stress tensor variables		 
         real(dp) :: exx, eyy, ezz, exz, eyz, exy
         real(dp) :: C
         real(dp) :: strxx, stryy, strzz, strxz, stryz, strxy
         real(dp) :: sxx, syy, szz, sxz, syz, sxy
       
         ! Per-voxel Full-Field Outputs for Post-Processing
         real(dp), allocatable  :: vm(:), stress_field(:,:)
       
         ! Microstructural state and iteration tracking
         real(dp), allocatable  :: orientation(:)      ! Euler angles per grain
         real(dp), allocatable  :: rotatedstiffness(:) ! Transformed stiffness tensor per grain
         real(dp), allocatable  :: energies(:), ggs(:) ! System energy tracking over iterations
         integer, allocatable  :: cgiter(:)            ! Iteration counter history

       end type elas3d_data_type
       
      end module elas3d_mod
    
!=========================================================================    
      program elas3d
      ! ==============================================================================
      ! PROGRAM: 3D Voxel-Based Linear Elasticity Solver
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Drives the finite element analysis to compute the effective macroscopic 
      !   elastic properties and local stress/strain fields of 3D microstructures 
      !   (e.g., composites, porous media, or polycrystals). 
      !
      ! MAIN PROGRAM VARIABLES (Conceptual Data Flow):
      !   - phase_map_file (intent(in)): 3D voxel microstructure input data.
      !   - material_moduli (intent(in)): Bulk/Shear moduli or anisotropic tensors.
      !   - nx, ny, nz (intent(in)): Grid dimensions defining the RVE size.
      !   - effective_stiffness (intent(out)): The homogenized macroscopic stiffness 
      !     tensor calculated at the end of the simulation.
      !   - local_fields (intent(out)): The resolved equilibrium displacement, stress, 
      !     and strain fields across the voxel grid.
      !
      ! ALGORITHM WORKFLOW:
      !   1. Initialization: Read inputs (microstructure, moduli) into the data struct.
      !   2. Pre-processing: Set up periodic neighbor mappings (ib array).
      !   3. FEM Assembly: Pass data to `femat_OpenMP` to generate solver inputs 
      !      (stiffness matrices, force vectors).
      !   4. Matrix Solve: Execute a Preconditioned Conjugate Gradient (PCG) solver 
      !      using the assembled outputs from step 3 to find the equilibrium 
      !      displacement field.
      !   5. Post-processing: Calculate local stresses/strains and homogenize to 
      !      find the effective properties.
      ! ==============================================================================	  
      use elas3d_mod
      use omp_lib   
      implicit none

      ! ======================================================================
      ! ---- 1. Variable Declarations ----------------------------------------
      ! ======================================================================       
      integer :: i1, j1, k1, m1
      integer :: i, j, k, m, n
      integer :: micro, npoints, ltot, kkk, Lstep
      integer :: ios, nxy
      integer :: im(nfaces), jm(nfaces), km(nfaces)
      real(dp) :: gg, utot, youngs
      real(dp) :: x, y, z, elapsed_time, gg0, phmod(2,2)
	  real(dp) :: rel_res
      real(dp), allocatable :: phasemod(:,:), prob(:)
      integer :: t1, t2, tc 
      integer :: lhist
      
      ! Instance of the global data structure defined in elas3d_mod	  
      type(elas3d_data_type) :: e3d
      
      
      ! ======================================================================
      ! ---- 2. Memory Allocation --------------------------------------------
      ! Allocates massive arrays dynamically based on parameters in module.
      ! ======================================================================    
      allocate(phasemod(nphase,2), prob(nphmax))
      allocate(e3d%u(ns,ndof), e3d%gb(ns,ndof))
      allocate(e3d%b(ns,ndof), e3d%h(ns,ndof))

      
      allocate(e3d%cmod(nphmax,6,6))
      allocate(e3d%dk(nphmax,nnode_fe,ndof,nnode_fe,ndof))
      allocate(e3d%vm(ns), e3d%stress_field(ns,6))

      allocate(e3d%ib(ns,nfaces), e3d%pix(ns))
      allocate(e3d%orientation(3*(nphase-1)))
      allocate(e3d%rotatedstiffness(36*(nphase-1)))    
      allocate(e3d%energies(n_iter+1))
      allocate(e3d%ggs(n_iter+1))
      allocate(e3d%cgiter(n_iter+1))

      ! ======================================================================
      ! ---- 3. Variable Initialization --------------------------------------
      ! Zeros out strain/stress accumulators and history arrays.
      ! ====================================================================== 
      e3d%exx = 0.d0;      e3d%eyy = 0.d0;      e3d%ezz = 0.d0;
      e3d%exz = 0.d0;      e3d%eyz = 0.d0;      e3d%exy = 0.d0;
      e3d%C   = 0.d0
      e3d%strxx = 0.d0;    e3d%stryy = 0.d0;    e3d%strzz = 0.d0;
      e3d%strxz = 0.d0;    e3d%stryz = 0.d0;    e3d%strxy = 0.d0;
      e3d%sxx = 0.d0;      e3d%syy = 0.d0;      e3d%szz = 0.d0;
      e3d%sxz = 0.d0;      e3d%syz = 0.d0;      e3d%sxy = 0.d0;
          
      ! CG iteration histories
      e3d%energies = 0.d0
      e3d%ggs      = 0.d0
      e3d%cgiter   = 0
      
      ! ======================================================================
      ! ---- 4. Setup Output Files and Start Clock ---------------------------
      ! ====================================================================== 
      open(unit=7, file='output_SS316L_polycrystal.txt', status='replace', iostat=ios)
      
      ! Timing Start (Fortran: use `system_clock`)
      call system_clock(t1, tc)  
      
      ! Output grid info to file and console
      write(7, '(a,i4,a,i4,a,i4,a,i12)') 'nx = ', nx, ', ny =', ny, ', nz =', nz, ', ns =', ns
      print '(a,i4,a,i4,a,i4,a,i12)', 'nx = ', nx, ', ny =', ny, ', nz =', nz, ', ns =', ns


      ! ======================================================================
      ! ---- 5. Initialize Material Properties (Isotropic Baseline) ----------
      ! Sets up the bulk and shear moduli for the matrix and inclusions.
      ! If flag_m=0, all matrix phases use Voigt-Reuss-Hill averages.
      ! If flag_m=1, only sets the inclusion (void) properties here.
      ! ====================================================================== 
      ! phasemod(i,1) = Bulk Modulus (E temporarily stored here), phasemod(i,2) = Poisson's ratio
      if (flag_m == 0) then
        do i = 1, nphase-1
           phasemod(i,1) = E0   ! matrix Young's Modulus
           phasemod(i,2) = nu0  ! matrix Poisson's Ratio
        end do
      else
	    phasemod(nphase,1) = E1  ! Inclusion Young's Modulus
		phasemod(nphase,2) = nu1 ! Inclusion Poisson's Ratio
      end if
      
      ! Convert Young's Modulus (E) and Poisson's Ratio (nu) 
      ! into Bulk Modulus (K) and Shear Modulus (G)
	  if (flag_m == 0) then
        do i = 1, nphase-1
          youngs = phasemod(i,1)
          phasemod(i,1) = phasemod(i,1) / 3.d0 / (1.d0 - 2.d0*phasemod(i,2)) ! K = E / 3(1-2v)
          phasemod(i,2) = youngs / 2.d0 / (1.d0 + phasemod(i,2))             ! G = E / 2(1+v) 
        end do
      else
	    youngs = phasemod(nphase,1)
		phasemod(nphase,1) = phasemod(nphase,1) / 3.d0 / (1.d0 - 2.d0*phasemod(nphase,2))
		phasemod(nphase,2) = youngs / 2.d0 / (1.d0 + phasemod(nphase,2))       
      end if
		     
      ! ======================================================================
      ! ---- 6. Build Periodic Boundary Neighbor Table -----------------------
      ! Constructs a 27-node local neighborhood map for every single node.
      ! Uses modulo arithmetic to wrap boundaries for periodic conditions.
      ! Optimized with OpenMP parallelization.
      ! ====================================================================== 
      ! Define relative 3D coordinate shifts for 27 neighbors
      im(1:8) = [0, 1, 1, 1, 0, -1, -1, -1]
      jm(1:8) = [1, 1, 0,-1,-1, -1,  0,  1]

      do n = 1, 8
        km(n)    = 0
        km(n+8)  = -1
        km(n+16) = 1
        im(n+8)  = im(n)
        im(n+16) = im(n)
        jm(n+8)  = jm(n)
        jm(n+16) = jm(n)
      end do

      im(25:27) = 0
      jm(25:27) = 0
      km(25) = -1; km(26) = 1; km(27) = 0      

      ! Map relative shifts to global 1-D array indices (with periodic wrapping)
      nxy = nx * ny
      !$omp parallel do collapse(3) schedule(guided) default(shared) &
      !$omp private(i,j,k,n,m,m1,i1,j1,k1)	  
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            m = nxy * (k-1) + nx * (j-1) + i
            do n = 1, 27
              i1 = i + im(n)
              j1 = j + jm(n)
              k1 = k + km(n)  
              ! Apply periodic wrapping			  
              if (i1 < 1)  i1 = i1 + nx
              if (i1 > nx) i1 = i1 - nx
              if (j1 < 1)  j1 = j1 + ny
              if (j1 > ny) j1 = j1 - ny
              if (k1 < 1)  k1 = k1 + nz
              if (k1 > nz) k1 = k1 - nz
              ! Compute 1D index of the wrapped neighbor			  
              m1 = nxy * (k1-1) + nx*(j1-1) + i1
              e3d%ib(m,n) = m1
            end do
          end do
        end do
      end do
      !$omp end parallel do
      
      ! ======================================================================
      ! ---- 7. Main Pre-processing Loop -------------------------------------
      ! Loads microstructure, calculates volume fractions, sets macroscopic 
      ! loads, builds stiffness matrices, and sets initial displacements.
      ! ====================================================================== 
      npoints = 1  ! Number of distinct microstructures to process
      do micro = 1, npoints
	  
        ! Read phase mapping and orientations from HDF5 file 
        call ppixel_hdf5(e3d)          
      
        print *, 'Unique phases in pix:', minval(e3d%pix), maxval(e3d%pix)
        do i = 1, nphase
          write(7, '(A,I7,A,F20.6,A,F20.6)') 'Phase ',i,' bulk=',phasemod(i,1),' shear=',phasemod(i,2)
        end do        
        
        ! Count and output the volume fractions of the different phases
		call assig(e3d, prob)
        do i = 1, nphase
          write(7, '(A,I7,A,F20.6)') 'Volume fraction of phase ',i,' is ',prob(i)
        end do

        ! Set applied macroscopic strains (defined in elas3d_mod)
        e3d%exx = aml_exx
        e3d%eyy = aml_eyy
        e3d%ezz = aml_ezz
        e3d%exz = aml_exz
        e3d%eyz = aml_eyz
        e3d%exy = aml_exy   
     
        write(7, *) 'Applied engineering strains'
        write(7, *) ' exx eyy ezz exz eyz exy'
        write(7, *) e3d%exx, e3d%eyy, e3d%ezz, 2.d0*e3d%exz, 2.d0*e3d%eyz, 2.d0*e3d%exy
        
        print *, 'Applied engineering strains: exx eyy ezz exz eyz exy'
        print '(6E15.5)', e3d%exx, e3d%eyy, e3d%ezz, 2.d0*e3d%exz, 2.d0*e3d%eyz, 2.d0*e3d%exy     


        ! Setup elastic modulus matrices, FE stiffness matrices, 
        ! constant C, and force vector b required for energy computation.  
		if (flag_m == 0) then
		  phmod(1:2,1:2) = phasemod(1:2,1:2)
		else
          phmod(1:1,1:2) = phasemod(nphase:nphase,1:2)    
        end if 
        ! Core matrix assembly routine (OpenMP parallelized)		
        call femat_OpenMP(e3d, phmod)
		
        ! Free up memory once stiffness is assembled		
        deallocate(e3d%orientation, e3d%rotatedstiffness)
      
        ! Apply macroscopic strains to create an initial guess for nodal displacements.   
        !$omp parallel do collapse(3) schedule(guided) default(shared) private(i,j,k,m,x,y,z)		
        do k = 1, nz
		  do j = 1, ny
            do i = 1, nx
			  m = nxy*(k-1)+nx*(j-1)+i
			  x = dble(i-1)
			  y = dble(j-1)
			  z = dble(k-1)
		   
              e3d%u(m,1) = x * e3d%exx + y * e3d%exy + z * e3d%exz
              e3d%u(m,2) = x * e3d%exy + y * e3d%eyy + z * e3d%eyz
              e3d%u(m,3) = x * e3d%exz + y * e3d%eyz + z * e3d%ezz
            end do
          end do
        end do
		!$omp end parallel do
     
        
      ! ======================================================================
      ! ---- 8. PCG Energy Minimization (Relaxation) Loop --------------------
      ! Iteratively solves the system Kd = F by minimizing the energy residual.
      ! ====================================================================== 
        ltot = 0
        lhist = 1
        
        ! Calculate initial energy and gradient using the initial displacement guess
        call energy_OpenMP(e3d, utot)
        
        ! Calculate initial gradient norm squared (used for relative tolerance check)
        gg = sum(e3d%gb(:,:)**2)
        gg0 = gg    ! Store initial norm as baseline
		
        write(7,*) 'Initial energy=', utot, 'gg=', gg
        write(7,*) '   '
        call flush(7)
        print '(A,E12.5,A,E12.5)', 'Initial energy=', utot, ', Initial gg=', gg
        print *, '   '
        
        e3d%energies(lhist) = utot
        e3d%ggs(lhist)      = sqrt(gg/gg0)
        e3d%cgiter(lhist)   = 0
        
        ! Open history file to track solver convergence
        open(unit=99, file='cg_history_polycrystal.txt', status='replace')
        write(99, '(I8,1X,1PE18.10,1X,1PE18.10)') e3d%cgiter(lhist), e3d%energies(lhist), e3d%ggs(lhist)
        
        ! Outer solver restart loop		
        do kkk = 1, kmax
		
           ! Core Conjugate Gradient Solver routine
           call dembx_OpenMP(e3d, Lstep, gg, gg0, kkk) 
                       
           ltot = ltot + Lstep
           lhist = lhist + 1
           
           ! Compute updated system energy after CG steps 
           call energy_OpenMP(e3d, utot)
         

           write(7,*) 'Energy=', utot, 'gg=', gg
           write(7,*) 'Number of conjugate steps=', ltot
           call flush(7) 
           print '(A,E12.5,A,E12.5,A,E12.5)', 'Energy=', utot, ', gg=', gg, ', rel. residual norm=', sqrt(gg/gg0)
           print *, 'Number of conjugate steps=', ltot
           
           ! Save CG history
           e3d%energies(lhist) = utot
           e3d%ggs(lhist)      = sqrt(gg/gg0)
           e3d%cgiter(lhist)   = ltot         
           write(99, '(I8,1X,1PE18.10,1X,1PE18.10)') e3d%cgiter(lhist), e3d%energies(lhist), e3d%ggs(lhist)

           ! Convergence Check: Jump out of loop if tolerance is met
		   rel_res = sqrt(gg/gg0)
           if (rel_res < tol) exit
           
           ! Compute and output intermediate macroscopic stresses/strains   
           call stress_OpenMP(e3d)
           
           write(7,*) 'stresses: xx, yy, zz, xz, yz, xy'
           write(7,'(6E15.6)') e3d%strxx, e3d%stryy, e3d%strzz, e3d%strxz, e3d%stryz, e3d%strxy
           write(7,*) 'strains: xx, yy, zz, xz, yz, xy'
           write(7,'(6E15.6)') e3d%sxx, e3d%syy, e3d%szz, e3d%sxz, e3d%syz, e3d%sxy
           write(7,*) '     '
           call flush(7)
           
        end do
        
        ! Close file for CG history
        close(99)

        ! Final computation of macroscopic stresses/strains
        call stress_OpenMP(e3d)

        write(7,*) 'stresses: xx, yy, zz, xz, yz, xy'
        write(7,'(6E15.6)') e3d%strxx, e3d%stryy, e3d%strzz, e3d%strxz, e3d%stryz, e3d%strxy
        write(7,*) 'strains: xx, yy, zz, xz, yz, xy'
        write(7,'(6E15.6)') e3d%sxx, e3d%syy, e3d%szz, e3d%sxz, e3d%syz, e3d%sxy

        
      end do
      
      print *, 'End of CG'
      
      ! ======================================================================
      ! ---- 9. Post-Processing and Teardown ---------------------------------
      ! Stops clock, calculates full-field stresses, and frees all memory.
      ! ====================================================================== 
      call system_clock(t2)
      elapsed_time = dble(t2-t1)/dble(tc)
      print *, ''
      print *, 'Calculation complete!'
      print *, 'Total calculation time: ', elapsed_time, ' seconds'   
	  write(7,*) ''
	  write(7,*) 'Calculation complete!'
      write(7,*) 'Total calculation time: ', elapsed_time, ' seconds'   	  
      
      ! Calculate and output voxel-by-voxel stress tensors and von Mises
      call stress_fullfield_OpenMP(e3d)    

      ! Free all dynamically allocated memory
      deallocate(phasemod, prob)
      deallocate(e3d%u, e3d%gb, e3d%b)
      deallocate(e3d%h)
      deallocate(e3d%cmod, e3d%dk)
      deallocate(e3d%vm, e3d%stress_field )
      deallocate(e3d%ib, e3d%pix)
      deallocate(e3d%energies, e3d%ggs, e3d%cgiter)

      
      end program elas3d    
    
!=========================================================================
      subroutine femat_OpenMP(e3d, phasemod)    
      ! ==============================================================================
      ! SUBROUTINE: femat_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Constructs the global finite element arrays required for the PCG solver. 
      !   It computes the local stiffness matrix for each 8-node hexahedral voxel 
      !   and applies periodic boundary conditions based on the prescribed macroscopic 
      !   strain state.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(inout)): The custom derived type containing the global FEM 
      !     arrays, grid data, and macroscopic boundary conditions.
      !   - phasemod (intent(in)): A 2x2 array holding the Bulk (K) and Shear (G) 
      !     moduli for the isotropic phases.
      !
      ! 'e3d' COMPONENTS MODIFIED/ACCESSED IN THIS SUBROUTINE:
      !   - e3d%pix (intent(in)): Voxel phase map indicating the material at each element.
      !   - e3d%ib (intent(in)): Periodic neighbor connectivity mapping.
      !   - e3d%exx, eyy, ezz... (intent(in)): Applied macroscopic strain components.
      !   - e3d%rotatedstiffness (intent(in)): Pre-calculated rotated tensors (polycrystals).
      !   - e3d%cmod (intent(out)): Populated with 6x6 local stiffness matrices.
      !   - e3d%dk (intent(out)): Populated with integrated finite element stiffness matrices.
      !   - e3d%b (intent(out)): Populated with the global right-hand-side force vector.
      !   - e3d%C (intent(out)): Populated with the global energy offset scalar.
      !
      ! METHODOLOGY:
      !   - Elements: 8-node trilinear hexahedral elements (voxels).
      !   - Integration: 3x3x3 Gauss quadrature (Simpson's rule integration weights).
      !   - Boundary Conditions: Enforces strict periodicity. 
      !   - Parallelization: Fully parallelized using OpenMP with atomic updates to 
      !     safely assemble the global force vector across multiple threads.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none

      ! ======================================================================
      ! ---- 1. Subroutine Arguments & Variable Declarations -----------------
      ! ====================================================================== 
      ! e3d: The global data structure containing state variables
      ! phasemod: 2x2 array holding Bulk (1) and Shear (2) moduli. 
      !           If polycrystal, it only holds the inclusion properties.
      type(elas3d_data_type), intent(inout) :: e3d
      real(dp), intent(in) :: phasemod(2,2)

      ! Local Loop Counters and Indexing Variables
      integer :: i, j, k, l, m, n
      integer :: ijk, mm, nn, ii, jj, kk, ll
      integer :: nxy, i3, i8, j3, k3, n8, m3, m8
	  
      ! Local Arrays for Element Formulation	  
      integer :: is(nnode_fe), ind
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
      real(dp) :: ck(6,6), cmu(6,6), g(ngauss,ngauss,ngauss)
      real(dp) :: es(6,nnode_fe,ndof), delta(nnode_fe,ndof)
      real(dp) :: left(6), right(6)
      real(dp) :: sumval, C
      real(dp) :: x, y, z
      real(dp) :: exx, eyy, ezz, exz, eyz, exy
	  
      ! Dynamically allocated local arrays for global assembly	  
      real(dp), allocatable :: dk(:,:,:,:,:)
      real(dp), allocatable :: b(:,:)

      ! ======================================================================
      ! ---- 2. Memory Allocation & Initialization ---------------------------
      ! Allocate thread-shared temporary arrays and strictly zero out all 
      ! arrays to prevent garbage values from contaminating the stiffness matrix.
      ! ====================================================================== 
      allocate(dk(nphmax,nnode_fe,ndof,nnode_fe,ndof))
	  allocate(b(ns,ndof))

      g = 0.d0;   es = 0.d0;   delta = 0.d0
      dk = 0.d0;  ck = 0.d0;   cmu = 0.d0
      left = 0.d0; right = 0.d0; b = 0.d0; C = 0.d0
      dndx = 0.d0; dndy = 0.d0; dndz = 0.d0
      e3d%cmod = 0.d0;   e3d%dk = 0.d0;   e3d%b = 0.d0
      e3d%C = 0.d0

      ! ======================================================================
      ! ---- 3. Finite Element Node Mapping ----------------------------------
      ! Maps the 8 corners of the hexahedral voxel element to specific indices 
      ! within the 27-node periodic neighbor array (ib) defined in main.
      ! ====================================================================== 
      is = (/27, 3, 2, 1, 26, 19, 18, 17/)	  
      nxy = nx * ny


      ! ======================================================================
      ! ---- 4. Isotropic Projection Matrices --------------------------------
      ! Sets up the volumetric (ck) and deviatoric (cmu) matrices used to 
      ! construct the 6x6 isotropic stiffness matrix from K and G moduli.
      ! ck matrix: Represents the bulk modulus (K) part.
      !            It’s an all-ones 3x3 matrix in the upper-left, zeros elsewhere.
      !            It contributes to normal stress/strain couplings.
      ! cmu matrix: Represents the shear (G) part.
      !             The upper-left 3x3 block (cmu(i,j) for i,j=1:3): 
      !             4/3 on diagonal, -2/3 off-diagonal.
      !             The lower-right diagonal blocks (cmu(4,4), cmu(5,5), cmu(6,6)): 1.0   
      ! ======================================================================
      ck(1:3,1:3) = 1.d0
      do i = 1, 3
        cmu(i,i) = 4.d0/3.d0
        do j = 1, 3
          if (i /= j) cmu(i,j) = -2.d0/3.d0
        end do
      end do
      cmu(4,4) = 1.d0
      cmu(5,5) = 1.d0
      cmu(6,6) = 1.d0      

      ! ======================================================================
      ! ---- 5. Constitutive Model Assembly (cmod) ---------------------------
      ! Builds the 6x6 local stiffness tensor for every phase in the domain.
      ! Compose cmod for each phase
      ! In Voigt notation for an isotropic material, the stiffness tensor can be constructed 
      ! using bulk modulus K, and shear modulus G, as follows:
      ! C_{ijkl} = K * del_{ij} * del_{kl} + 
      !            2G(1/2*(del_{ik}*del_{jl}+del_{il}*del_{jk}) - 1/3*(del_{ij}*del_{kl}))
      ! In 6x6 Voigt matrix form, this can be written as:
      ! C_{ij} = K * ck_{ij} + G * cmu_{ij} 
      ! ====================================================================== 	  
	  if (flag_m == 0) then
        ! Case A: Isotropic Matrix 
        ! Combines volumetric and deviatoric parts for phases 1 & 2	  
	    do k = 1,2; do j = 1,6; do i = 1,6;
		  e3d%cmod(k,i,j) = phasemod(k,1)*ck(i,j) + phasemod(k,2)*cmu(i,j)
        end do; end do; end do		
      else
        ! Case B: Polycrystalline Matrix
        ! Computes Euler rotations to transform single-crystal stiffness to global axes	  
		call rotatestiffnessby(e3d)	
		
        ! Assign the transformed anisotropic stiffness to each grain phase
        !$omp parallel do collapse(3) schedule(guided) default(shared)
        do k = 1, nphase-1
		  do i = 1,6
		    do j = 1,6
              ! Maps the 1D flattened rotated stiffness array into a 3D array			
			  e3d%cmod(k,i,j) = e3d%rotatedstiffness(36*(k-1)+(i-1)*6+j)
			end do
          end do
        end do
		!$omp end parallel do
      
        ! Assign isotropic properties (0.0 stiffness) to the inclusion/pore phase
        do j = 1, 6
		  do i = 1, 6
            e3d%cmod(nphase,i,j) = phasemod(1,1)*ck(i,j) + phasemod(1,2)*cmu(i,j)
          end do
		end do   
      end if

      ! ======================================================================
      ! ---- 6. Numerical Integration Weights (Simpson's Rule) ---------------
      ! Sets up a 3x3x3 grid of weighting factors (1, 4, 16, 64) for numerical 
      ! integration of the finite element stiffness matrices.
      ! ====================================================================== 
      do k3 = 1, ngauss; do j3 = 1, ngauss; do i3 = 1, ngauss;
        n = 0
        ! Count how many 'middle' nodes this point corresponds to		
        if(i3==2) n = n+1
        if(j3==2) n = n+1
        if(k3==2) n = n+1
        ! Assign weight based on dimensionality (4^n)		
        g(i3,j3,k3) = 4.0 ** n
      end do; end do; end do

      ! ======================================================================
      ! ---- 7. Main Element Stiffness Assembly (Gauss Quadrature) -----------
      ! Integrates the stiffness matrix (K = B^T * D * B) over the volume of
      ! the hexahedral element for each phase using 3x3x3 Gauss integration.
      ! ======================================================================
      !$omp parallel default(shared) private(k3,j3,i3,x,y,z,dndx,dndy,dndz,es,n8,mm,nn,ii,jj,left,right,sumval)
      !$omp do schedule(guided)
      do ijk = 1, nphase
	    do k3 = 1, ngauss
          do j3 = 1, ngauss
            do i3 = 1, ngauss;
              ! Local element coordinates (x, y, z) mapped from 0 to 1			
              x = dble(i3-1)/2.d0
              y = dble(j3-1)/2.d0
              z = dble(k3-1)/2.d0

              ! Derivatives of the 8-node trilinear shape functions (dN/dx, dN/dy, dN/dz)
              dndx(1) = -(1.d0-y)*(1.d0-z)
              dndx(2) =  (1.d0-y)*(1.d0-z)
              dndx(3) =         y*(1.d0-z)
              dndx(4) =        -y*(1.d0-z)
              dndx(5) = -(1.d0-y)*z
              dndx(6) =  (1.d0-y)*z
              dndx(7) =         y*z
              dndx(8) =        -y*z

              dndy(1) = -(1.d0-x)*(1.d0-z)
              dndy(2) =        -x*(1.d0-z)
              dndy(3) =         x*(1.d0-z)
              dndy(4) =  (1.d0-x)*(1.d0-z)
              dndy(5) = -(1.d0-x)*z
              dndy(6) =        -x*z
              dndy(7) =         x*z
              dndy(8) =  (1.d0-x)*z

              dndz(1) = -(1.d0-x)*(1.d0-y)
              dndz(2) =        -x*(1.d0-y)
              dndz(3) =        -x*y
              dndz(4) = -(1.d0-x)*y
              dndz(5) =  (1.d0-x)*(1.d0-y)
              dndz(6) =         x*(1.d0-y)
              dndz(7) =         x*y
              dndz(8) =  (1.d0-x)*y

              ! Construct the Strain-Displacement Matrix (B-matrix), called `es` here
              ! Maps nodal degrees of freedom (u, v, w) to 6 strain components
              es = 0.d0
              do n8 = 1, nnode_fe
                es(1,n8,1) = dndx(n8)
                es(2,n8,2) = dndy(n8)
                es(3,n8,3) = dndz(n8)
                es(4,n8,1) = dndz(n8)
                es(4,n8,3) = dndx(n8)
                es(5,n8,2) = dndz(n8)
                es(5,n8,3) = dndy(n8)
                es(6,n8,1) = dndy(n8)
                es(6,n8,2) = dndx(n8)
              end do

              ! Perform matrix multiplication: B^T * D * B
              ! Accumulate into local stiffness matrix `dk` using Gauss weights
              do mm = 1, ndof
                do nn = 1, ndof
                  do ii = 1, nnode_fe
                    left = es(:, ii, mm)
                    do jj = 1, nnode_fe
                      right = es(:, jj, nn)
                      sumval = dot_product(left, matmul(e3d%cmod(ijk,:,:), right))
                      ! Divide by 216 to normalize the integration volume					  
                      dk(ijk,ii,mm,jj,nn) = dk(ijk,ii,mm,jj,nn) + g(i3,j3,k3)*sumval/216.d0
                    end do
                  end do
                end do
              end do

            end do
          end do
        end do
      end do
	  !$omp end do
	  !$omp end parallel

      e3d%dk = dk ! Assign thread-local completed matrix to global struct

      ! ======================================================================
      ! ---- 8. Periodic Boundary Conditions & Force Vectors -----------------
      ! To enforce a uniform macroscopic strain on a periodic grid, nodes on
      ! the upper boundary faces, edges, and corners must have a displacement
      ! "jump" (delta) applied relative to their periodic counterparts.
      ! This creates a body force equivalent (b) and energy offset (C).
      ! ======================================================================
	  
      b = e3d%b
      C = e3d%C
	  
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy

      ! -------- Face: x = nx --------
      ! Calculate displacement jump (delta) for nodes on this face
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==2 .or. i8==3 .or. i8==6 .or. i8==7) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
      end do

      ! Apply force corrections due to displacement jump
      !$omp parallel do collapse(2) schedule(guided) private(j,k,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do j = 1, ny-1
        do k = 1, nz-1
          m = nxy*(k-1) + j*nx
          do nn = 1, ndof
            do mm = 1, nnode_fe
              sumval = 0.d0
              do m3 = 1, ndof
                do m8 = 1, nnode_fe
                  sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
                  C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn) * delta(mm,nn)
                end do
              end do
              ind = e3d%ib(m,is(mm))			  
              ! Use atomic update to avoid thread collisions on the global b vector			  
              !$omp atomic update
			  b(ind,nn) = b(ind,nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Face: y = ny --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==3 .or. i8==4 .or. i8==7 .or. i8==8) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
      end do;

      !$omp parallel do collapse(2) schedule(guided) private(i,k,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do i = 1, nx-1
        do k = 1, nz-1
          m = nxy*(k-1) + nx*(ny-1) + i
          do nn = 1, ndof
            do mm = 1, nnode_fe
              sumval = 0.d0
              do m3 = 1, ndof
                do m8 = 1, nnode_fe
                  sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
                  C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn) * delta(mm,nn)
                end do
              end do
			  ind = e3d%ib(m,is(mm))
              !$omp atomic update
			  b(ind,nn) = b(ind,nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Face: z = nz --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==5 .or. i8==6 .or. i8==7 .or. i8==8) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
      end do;

      !$omp parallel do collapse(2) schedule(guided) private(i,j,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do i = 1, nx-1
        do j = 1, ny-1
          m = nxy*(nz-1) + nx*(j-1) + i
          do nn = 1, ndof
            do mm = 1, nnode_fe
              sumval = 0.d0
              do m3 = 1, ndof
                do m8 = 1, nnode_fe
                  sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
                  C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn) * delta(mm,nn)
                end do
              end do
              ind = e3d%ib(m,is(mm))
              !$omp atomic update
			  b(ind,nn) = b(ind,nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Edge: x=nx, y=ny --------
      ! Nodes on edges require combined displacement jumps from both adjacent faces
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==2 .or. i8==6) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
        if (i8==4 .or. i8==8) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
        if (i8==3 .or. i8==7) then
          delta(i8,1) = exy*ny + exx*nx
          delta(i8,2) = eyy*ny + exy*nx
          delta(i8,3) = eyz*ny + exz*nx
        end if
      end do;

      !$omp parallel do schedule(guided) private(k,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do k = 1, nz-1
        m = nxy*k
        do nn = 1, ndof
          do mm = 1, nnode_fe
            sumval = 0.d0
            do m3 = 1, ndof
              do m8 = 1, nnode_fe
                sumval = sumval + delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)
                C = C + 0.5d0 * delta(m8,m3) * dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
              end do
            end do
            ind = e3d%ib(m,is(mm))
            !$omp atomic update
            b(ind,nn) = b(ind,nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Edge: x=nx, z=nz --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if(i8==2 .or. i8==3) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
        if(i8==5 .or. i8==8) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
        if(i8==6 .or. i8==7) then
          delta(i8,1) = exz*nz + exx*nx
          delta(i8,2) = eyz*nz + exy*nx
          delta(i8,3) = ezz*nz + exz*nx
        end if
      end do;

      !$omp parallel do schedule(guided) private(j,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do j = 1, ny-1
        m = nxy*(nz-1) + nx*(j-1) + nx
        do nn = 1, ndof
          do mm = 1, nnode_fe
            sumval = 0.d0
            do m3 = 1, ndof
              do m8 = 1, nnode_fe
                sumval = sumval + delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)
                C = C + 0.5d0 * delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
              end do
            end do
            ind = e3d%ib(m,is(mm))
            !$omp atomic update
            b(ind,nn) = b(ind,nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Edge: y=ny, z=nz --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if(i8==5 .or. i8==6) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
        if(i8==3 .or. i8==4) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
        if(i8==7 .or. i8==8) then
          delta(i8,1) = exy*ny + exz*nz
          delta(i8,2) = eyy*ny + eyz*nz
          delta(i8,3) = eyz*ny + ezz*nz
        end if
      end do;

      !$omp parallel do schedule(guided) private(i,m,nn,mm,sumval,m3,m8) reduction(+:C)	  
      do i = 1, nx-1
        m = nxy*(nz-1) + nx*(ny-1) + i
        do nn = 1, ndof
          do mm = 1, nnode_fe
            sumval = 0.d0
            do m3 = 1, ndof
              do m8 = 1, nnode_fe
                sumval = sumval + delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)
                C = C + 0.5d0 * delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
              end do
            end do
            ind = e3d%ib(m,is(mm))
            !$omp atomic update
            b(ind,nn) = b(ind,nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- Corner: x=nx, y=ny, z=nz --------
      ! The single corner node requires cumulative jumps from all 3 dimensions
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if(i8==2) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
        if(i8==4) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
        if(i8==5) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
        if(i8==8) then
          delta(i8,1) = exy*ny + exz*nz
          delta(i8,2) = eyy*ny + eyz*nz
          delta(i8,3) = eyz*ny + ezz*nz
        end if
        if(i8==6) then
          delta(i8,1) = exx*nx + exz*nz
          delta(i8,2) = exy*nx + eyz*nz
          delta(i8,3) = exz*nx + ezz*nz
        end if
        if(i8==3) then
          delta(i8,1) = exx*nx + exy*ny
          delta(i8,2) = exy*nx + eyy*ny
          delta(i8,3) = exz*nx + eyz*ny
        end if
        if(i8==7) then
          delta(i8,1) = exx*nx + exy*ny + exz*nz
          delta(i8,2) = exy*nx + eyy*ny + eyz*nz
          delta(i8,3) = exz*nx + eyz*ny + ezz*nz
        end if
      end do;

      m = nx*ny*nz
      do nn = 1, ndof
        do mm = 1, nnode_fe
          sumval = 0.d0
          do m3 = 1, ndof
            do m8 = 1, nnode_fe
              sumval = sumval + delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)
              C = C + 0.5d0 * delta(m8,m3)*dk(e3d%pix(m),m8,m3,mm,nn)*delta(mm,nn)
            end do
          end do
          b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
        end do
      end do
	  
      ! ======================================================================
      ! ---- 9. Final Variable Assignment & Cleanup --------------------------
      ! ======================================================================
      e3d%b = b
      e3d%C = C
   
      deallocate(dk, b)

      end subroutine femat_OpenMP

!-------------------------------------------------------------------------   
      subroutine energy_OpenMP(e3d, utot)
      ! ==============================================================================
      ! SUBROUTINE: energy_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Computes the total strain energy of the system and calculates the energy 
      !   gradient vector. This subroutine acts as the matrix-free matrix-vector 
      !   multiplication step for the iterative solver (e.g., Conjugate Gradient), 
      !   implicitly calculating K * u without assembling the global stiffness matrix K.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(inout)): The custom derived type containing the global FEM arrays.
      !   - utot (intent(out)): Total calculated strain energy of the current state.
      !
      ! 'e3d' COMPONENTS ACCESSED/MODIFIED (Explicit Intent Meaning):
      !   - e3d%C (intent(in)): Base energy constant offset from macroscopic boundaries.
      !   - e3d%u (intent(in)): Current equilibrium displacement field (ns x ndof).
      !   - e3d%b (intent(in)): Global right-hand-side equivalent force vector.
      !   - e3d%dk (intent(in)): Local integrated finite element stiffness matrices.
      !   - e3d%pix (intent(in)): 1D voxel phase map.
      !   - e3d%ib (intent(in)): 27-node periodic neighborhood connectivity mapping.
      !   - e3d%gb (intent(out)): Computed energy gradient vector (residual force).
      !
      ! METHODOLOGY:
      !   - Calculates internal forces by summing contributions from the 27 surrounding 
      !     neighborhood nodes for each voxel (x, y, z degrees of freedom).
      !   - Total Energy Formula: U = 0.5 * u^T * K * u + b^T * u + C
      !   - Gradient Formula: g = K * u + b
      ! ==============================================================================	  
      use elas3d_mod
      implicit none

      ! ======================================================================
      ! ---- 1. Subroutine Arguments & Variable Declarations -----------------
      ! ====================================================================== 
      ! e3d: The global data structure containing mesh mapping, physical 
      !      properties, displacements, and local element stiffness matrices.
      ! utot: Total energy or macroscopic output parameter being calculated.
      type(elas3d_data_type), intent(inout) :: e3d
      real(dp), intent(out) :: utot

      ! Local Loop Counters and Iteration Variables
      integer :: m, m3, j
      real(dp) :: C, gbsum
	  
      ! Dynamically allocated local array for gradient/force storage	  
      real(dp), allocatable :: gb(:,:)

      ! ======================================================================
      ! ---- 2. Memory Allocation & Initialization ---------------------------
      ! Allocates the gradient/internal force vector (gb) for all spatial 
      ! nodes (ns) and degrees of freedom (ndof). 
      ! ====================================================================== 
      allocate( gb(ns,ndof) )

      ! Unpack the macroscopic background energy offset from struct
      C   = e3d%C

      ! Strictly zero out the array to prevent garbage values from contaminating
      ! the force accumulations in the subsequent parallel loops.
      gb(:,:) = 0.d0
	  
      ! ======================================================================
      ! ---- 3. Global Force / Gradient Vector Assembly (Matrix-Free) --------
      ! Computes gb = A * u using a matrix-free sparse multiplication approach.
      ! Outer loop iterates over the spatial degrees of freedom (j = 1 to ndof).
      ! Inner loop iterates over all voxels/nodes (m = 1 to ns) and is 
      ! parallelized via OpenMP for high-performance 3D grid evaluation.
      ! ====================================================================== 	  
      do j = 1, ndof
        !$omp parallel do private(m, gbsum) schedule(guided)		
        do m = 1, ns
          gbsum = 0.0d0

          ! ==================================================================
          ! ---- 27-Node Hexahedral Element Stencil Computations (n=1) ----
          ! Calculates the internal force/energy gradient acting on voxel 'm'.
          ! This block strictly handles the X-direction (n=1) displacement 
          ! contributions from the voxel's 27 periodic neighbors.
          !
          ! Anatomy of the calculation:
          ! e3d%ib(m, idx) : Retrieves global index of neighbor 'idx' (1-27).
          ! e3d%u(..., 1)  : Retrieves the X-displacement of that neighbor.
          ! e3d%dk(...)    : Local stiffness mapped by the neighbor's material
          !                  phase (pix), local node idx, and current DoF (j).
          ! ================================================================== 
          ! === n = 1 ===
          gbsum = gbsum                                 &
        + e3d%u(e3d%ib(m,1),1 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,7,1) ) &
        + e3d%u(e3d%ib(m,2),1 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,7,1) )                                            &
        + e3d%u(e3d%ib(m,3),1 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,7,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,6,1) ) &
        + e3d%u(e3d%ib(m,4),1 )*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,6,1) )                                             &
        + e3d%u(e3d%ib(m,5),1 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,6,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,5,1) ) &
        + e3d%u(e3d%ib(m,6),1 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,5,1) )                                              &
        + e3d%u(e3d%ib(m,7),1 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,5,1) ) &
        + e3d%u(e3d%ib(m,8),1 )*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,8,1) )                                              &
        + e3d%u(e3d%ib(m,9),1 )*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,3,1) )                                             &
        + e3d%u(e3d%ib(m,10),1)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,3,1) )                                               &
        + e3d%u(e3d%ib(m,11),1)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,2,1) )                                             &
        + e3d%u(e3d%ib(m,12),1)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,2,1) )                                               &
        + e3d%u(e3d%ib(m,13),1)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,2,1) )                                             &
        + e3d%u(e3d%ib(m,14),1)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,1,1) )                                               &
        + e3d%u(e3d%ib(m,15),1)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,4,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,1,1) )                                             &
        + e3d%u(e3d%ib(m,16),1)*( e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,4,1) )                                               &
        + e3d%u(e3d%ib(m,17),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,7,1) )                                              &
        + e3d%u(e3d%ib(m,18),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,7,1) )                                               &
        + e3d%u(e3d%ib(m,19),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,6,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,7,1) )                                              &
        + e3d%u(e3d%ib(m,20),1)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,6,1) )                                                &
        + e3d%u(e3d%ib(m,21),1)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,5,1) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,6,1) )                                               &
        + e3d%u(e3d%ib(m,22),1)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,5,1) )                                                &
        + e3d%u(e3d%ib(m,23),1)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,8,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,5,1) )                                              &
        + e3d%u(e3d%ib(m,24),1)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,8,1) )                                                &
        + e3d%u(e3d%ib(m,25),1)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,4,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,2,1) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,1,1) ) &
        + e3d%u(e3d%ib(m,26),1)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,7,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,8,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,5,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,6,1) ) &
        + e3d%u(e3d%ib(m,27),1)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,1,1) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,2,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,3,1) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,4,1) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,5,1) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,6,1) + &
                                      e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,7,1) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,8,1) )

          ! ==================================================================
          ! ---- Y-Direction (n=2) Displacement Contributions ----------------
          ! Accumulates the internal force/energy gradient acting on voxel 'm'
          ! based on the Y-direction (n=2) displacements of its 27 periodic 
          ! neighbors.
          !
          ! Note: e3d%u(..., 2) extracts the Y-displacement of the neighbor,
          ! and e3d%dk(..., 2) accesses the corresponding Y-direction stiffness
          ! matrix components. This block structurally sits between n=1 and n=3.
          ! ================================================================== 
          ! === n = 2 ===
          gbsum = gbsum                                 &
        + e3d%u(e3d%ib(m,1),2 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,7,2) ) &
        + e3d%u(e3d%ib(m,2),2 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,7,2) )                                            &
        + e3d%u(e3d%ib(m,3),2 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,7,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,6,2) ) &
        + e3d%u(e3d%ib(m,4),2 )*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,6,2) )                                             &
        + e3d%u(e3d%ib(m,5),2 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,6,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,5,2) ) &
        + e3d%u(e3d%ib(m,6),2 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,5,2) )                                              &
        + e3d%u(e3d%ib(m,7),2 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,5,2) ) &
        + e3d%u(e3d%ib(m,8),2 )*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,8,2) )                                              &
        + e3d%u(e3d%ib(m,9),2 )*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,3,2) )                                             &
        + e3d%u(e3d%ib(m,10),2)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,3,2) )                                               &
        + e3d%u(e3d%ib(m,11),2)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,2,2) )                                             &
        + e3d%u(e3d%ib(m,12),2)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,2,2) )                                               &
        + e3d%u(e3d%ib(m,13),2)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,2,2) )                                             &
        + e3d%u(e3d%ib(m,14),2)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,1,2) )                                               &
        + e3d%u(e3d%ib(m,15),2)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,4,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,1,2) )                                             &
        + e3d%u(e3d%ib(m,16),2)*( e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,4,2) )                                               &
        + e3d%u(e3d%ib(m,17),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,7,2) )                                              &
        + e3d%u(e3d%ib(m,18),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,7,2) )                                               &
        + e3d%u(e3d%ib(m,19),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,6,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,7,2) )                                              &
        + e3d%u(e3d%ib(m,20),2)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,6,2) )                                                &
        + e3d%u(e3d%ib(m,21),2)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,5,2) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,6,2) )                                               &
        + e3d%u(e3d%ib(m,22),2)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,5,2) )                                                &
        + e3d%u(e3d%ib(m,23),2)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,8,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,5,2) )                                              &
        + e3d%u(e3d%ib(m,24),2)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,8,2) )                                                &
        + e3d%u(e3d%ib(m,25),2)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,4,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,2,2) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,1,2) ) &
        + e3d%u(e3d%ib(m,26),2)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,7,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,8,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,5,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,6,2) ) &
        + e3d%u(e3d%ib(m,27),2)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,1,2) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,2,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,3,2) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,4,2) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,5,2) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,6,2) + &
                                      e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,7,2) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,8,2) )

          ! ==================================================================
          ! ---- Z-Direction (n=3) Displacement Contributions -------------
          ! Completes the internal force accumulation for voxel 'm' by adding 
          ! the structural pull caused by the Z-direction (n=3) displacements 
          ! of its 27 periodic neighbors.
          ! 
          ! Note: e3d%u(..., 3) extracts the Z-displacement of the neighbor,
          ! and e3d%dk(..., 3) accesses the corresponding Z-direction stiffness.
          ! (Assuming n=2 was calculated prior to this block).
          ! ================================================================== 
          ! === n = 3 ===
          gbsum = gbsum                                 &
        + e3d%u(e3d%ib(m,1),3 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,7,3) ) &
        + e3d%u(e3d%ib(m,2),3 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,7,3) )                                            &
        + e3d%u(e3d%ib(m,3),3 )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,7,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,6,3) ) &
        + e3d%u(e3d%ib(m,4),3 )*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,6,3) )                                             &
        + e3d%u(e3d%ib(m,5),3 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,6,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,5,3) ) &
        + e3d%u(e3d%ib(m,6),3 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,5,3) )                                              &
        + e3d%u(e3d%ib(m,7),3 )*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,5,3) ) &
        + e3d%u(e3d%ib(m,8),3 )*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,8,3) )                                              &
        + e3d%u(e3d%ib(m,9),3 )*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,3,3) )                                             &
        + e3d%u(e3d%ib(m,10),3)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,3,3) )                                               &
        + e3d%u(e3d%ib(m,11),3)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,2,3) )                                             &
        + e3d%u(e3d%ib(m,12),3)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,2,3) )                                               &
        + e3d%u(e3d%ib(m,13),3)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,2,3) )                                             &
        + e3d%u(e3d%ib(m,14),3)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,1,3) )                                               &
        + e3d%u(e3d%ib(m,15),3)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,4,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,1,3) )                                             &
        + e3d%u(e3d%ib(m,16),3)*( e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,4,3) )                                               &
        + e3d%u(e3d%ib(m,17),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,7,3) )                                              &
        + e3d%u(e3d%ib(m,18),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,7,3) )                                               &
        + e3d%u(e3d%ib(m,19),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,6,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,7,3) )                                              &
        + e3d%u(e3d%ib(m,20),3)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,6,3) )                                                &
        + e3d%u(e3d%ib(m,21),3)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,5,3) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,6,3) )                                               &
        + e3d%u(e3d%ib(m,22),3)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,5,3) )                                                &
        + e3d%u(e3d%ib(m,23),3)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,8,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,5,3) )                                              &
        + e3d%u(e3d%ib(m,24),3)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,8,3) )                                                &
        + e3d%u(e3d%ib(m,25),3)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,4,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,2,3) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,1,3) ) &
        + e3d%u(e3d%ib(m,26),3)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,7,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,8,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,5,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,6,3) ) &
        + e3d%u(e3d%ib(m,27),3)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,1,3) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,2,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,3,3) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,4,3) +  &
                                      e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,5,3) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,6,3) + &
                                      e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,7,3) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,8,3) )

          gb(m, j) = gbsum
        end do
        !$omp end parallel do
      end do

      ! ======================================================================
      ! ---- 4. Total Strain Energy & Final Gradient Update ------------------
      ! Calculates the total macroscopic system energy (utot) using the 
      ! standard quadratic form: U = C + 0.5*(u^T * A * u) + (b^T * u).
      ! Here, gb initially holds the (A * u) term from the previous loops.
      ! It also finalizes the gradient vector (gb) by adding the linear 
      ! force/offset term (b). The OpenMP reduction(+:utot) clause safely 
      ! accumulates the scalar energy across all threads.
      ! ======================================================================
      utot = C
      do m3 = 1, ndof
        !$omp parallel do schedule(guided) reduction(+:utot)		
        do m = 1, ns
          utot = utot + 0.5d0 * e3d%u(m, m3) * gb(m, m3) + e3d%b(m, m3) * e3d%u(m, m3)
          gb(m, m3) = gb(m, m3) + e3d%b(m, m3)
        end do
        !$omp end parallel do
      end do

      ! ======================================================================
      ! ---- 5. Global State Update & Memory Cleanup -------------------------
      ! Copies the fully assembled local gradient array back into the global 
      ! data structure (e3d) so it can be accessed by the main solver loop 
      ! Finally, deallocates local arrays to prevent memory leaks upon 
      ! exiting the subroutine.
      ! ======================================================================
      e3d%gb = gb

      deallocate( gb )

      end subroutine energy_OpenMP
  
!-------------------------------------------------------------------------
      subroutine dembx_OpenMP(e3d, Lstep, gg, gg0, kkk)
      ! ==============================================================================
      ! SUBROUTINE: dembx_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Iterative solver that computes the equilibrium displacement field of the 
      !   microstructure. It solves the large, sparse system of linear equations 
      !   using a Preconditioned Conjugate Gradient (PCG) algorithm to minimize 
      !   the total elastic strain energy.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(inout)): The custom derived type containing the global FEM 
      !     arrays, including displacements, gradients, and stiffness components.
      !   - Lstep (intent(out)): The total number of iterations required to reach 
      !     convergence in the current load step.
      !   - gg (intent(out)): The final relative residual norm upon exiting the loop.
      !   - gg0 (intent(in)): The initial residual norm, used to calculate the 
      !     relative convergence tolerance.
      !   - kkk (intent(in)): The current load step index. When kkk == 1, it triggers 
      !     the assembly and inversion of the Block Jacobi preconditioner.
      !
      ! 'e3d' COMPONENTS MODIFIED/ACCESSED IN THIS SUBROUTINE:
      !   - e3d%ib (intent(in)): Periodic neighbor connectivity mapping (27 nodes).
      !   - e3d%pix (intent(in)): Voxel phase map indicating the material phase.
      !   - e3d%dk (intent(in)): Local finite element stiffness matrix components.
      !   - e3d%u (intent(inout)): The 3D displacement field vector (updated iteratively).
      !   - e3d%gb (intent(inout)): The residual/gradient vector (updated iteratively).
      !   - e3d%h (intent(inout)): The conjugate search direction vector.
      !
      ! METHODOLOGY:
      !   - Algorithm: Preconditioned Conjugate Gradient (PCG) method.
      !   - Preconditioner: 3x3 Block Jacobi. Extracts the diagonal 3x3 stiffness 
      !     block for each node (summed from its 8 adjacent octants) and applies 
      !     its inverse to pre-condition the residual, drastically reducing iterations.
      !   - Sparse Matrix-Vector Multiply (SpMV): Calculates the structural response 
      !     (Ah = A*h) using a cache-friendly blocked loop over the 27-node neighborhood.
      !   - Divergence Control: Monitors the relative residual. If it increases by 
      !     more than an order of magnitude (10x), the search direction is reset to 
      !     steepest descent to recover orthogonality and prevent divergence.
      !   - Parallelization: Fully parallelized using OpenMP (`schedule(guided)`). 
      !     Matrix operations, inner products, and vector updates are distributed 
      !     across threads using local fast-memory arrays to avoid race conditions.
      !   - I/O: Logs the convergence history (iteration vs. residual) to 
      !     'cgitr_SS316L_polycrystal.txt'.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none

      ! ======================================================================
      ! ---- 1. Subroutine Arguments -----------------------------------------
      ! Iterative Preconditioned Conjugate Gradient (PCG) Solver.
      ! e3d:   Global data structure containing the FEA state.
      ! Lstep: Total number of iterations taken to reach convergence.
      ! gg:    Current residual norm (error squared) upon exit.
      ! gg0:   Initial or target residual norm for convergence criteria.
      ! kkk:   Maximum allowed iterations or current load step index.
      ! ======================================================================
      type(elas3d_data_type), intent(inout) :: e3d
      integer, intent(out) :: Lstep
      real(dp), intent(out) :: gg
      real(dp), intent(in)  :: gg0
      integer, intent(in)   :: kkk

      ! ======================================================================
      ! ---- 2. Local Array Declarations -------------------------------------
      ! Dynamically allocated arrays to hold the solver's state. 
      ! Ah:          Stores the sparse matrix-vector product $A \mathbf{h}$.
      ! M_block_inv: Inverse of the 3x3 Block Jacobi preconditioner matrix.
      ! z:           The preconditioned residual vector $\mathbf{z} = \mathbf{M}^{-1} \mathbf{g}$.
      ! history_*:   Convergence tracking arrays for I/O logging.
      ! ====================================================================== 
      ! --- Local Arrays (Saved to avoid re-allocation cost) ---
      real(dp), allocatable :: u(:,:), gb(:,:), h(:,:)
      real(dp), allocatable :: Ah(:,:)
      real(dp), allocatable :: M_block_inv(:,:,:) ! Precomputed inverse
      real(dp), allocatable :: z(:,:)             ! Preconditioned residual  
      integer, allocatable :: history_step(:)
      real(dp), allocatable :: history_res(:)
        
      ! ======================================================================
      ! ---- 3. Local Scalars & CG Variables ---------------------------------
      ! Standard scalars for PCG math: lambda (step size $\alpha$), gamma 
      ! (orthogonalization factor $\beta$), and hAh (curvature $\mathbf{h}^T A \mathbf{h}$).
      ! ======================================================================
      integer :: ns_local, L1, L2
      integer :: m, ijk, i, j, k
      integer :: ios                       		  
      real(dp) :: lambda, gamma, hAh, gglast, rel_res
      real(dp) :: rel_res_prev, rel_res_min
      integer :: reset_count		  

	  
      ! ======================================================================
      ! ---- 4. Thread-Private Variables -------------------------------------
      ! Variables that must remain private to each OpenMP thread to avoid 
      ! race conditions during the parallel assembly and matrix inversion.
      ! ======================================================================
      integer :: ibl(27), pl(27)
      real(dp) :: ahsum1, ahsum2, ahsum3 
	  real(dp) :: M_block(3,3)             ! Block Jacobi preconditioner
	  
      ! Open file to log the Conjugate Gradient iterations for SS316L polycrystal	  
	  open(unit=100, file='cgitr_SS316L_polycrystal.txt', status='replace', iostat=ios)

      ns_local = ns

      ! ======================================================================
      ! ---- 5. Memory Allocation --------------------------------------------
      ! Allocates solver arrays using the global number of spatial nodes (ns)
      ! and degrees of freedom (ndof).
      ! ======================================================================
      allocate( u(ns,ndof), gb(ns,ndof), h(ns,ndof) )
      allocate( Ah(ns,ndof) )
	  allocate( M_block_inv(ns,ndof,ndof) )
      allocate( z(ns,ndof) )
      allocate(history_step(ldemb))
      allocate(history_res(ldemb))	  

      ! ======================================================================
      ! ---- 6. Local State Unpacking ----------------------------------------
      ! Copies the global displacements (u), gradient/residual (gb), and 
      ! conjugate search direction (h) into fast local memory arrays.
      ! Note: This relies on an enclosing !$omp parallel region, or implies 
      ! threaded execution initiated by the caller.
      ! ======================================================================
	  !$omp do schedule(guided)
      do m = 1, ns
        u(m,:)  = e3d%u(m,:)
        gb(m,:) = e3d%gb(m,:)
        h(m,:)  = e3d%h(m,:)
      end do
      !$omp end do
	  
  
      ! ======================================================================
      ! ---- 7. Block Jacobi Preconditioner Assembly (First Step Only) -------
      ! If this is the first load step (kkk == 1), we build and invert a 
      ! 3x3 preconditioner matrix (M_block) for every node. This block 
      ! represents the local diagonal of the global stiffness matrix A.
      ! It is assembled by summing the stiffness contributions from the 8 
      ! surrounding elements (octants) connected to node 'm'.
      ! ======================================================================
	  if (kkk == 1) then
	    ! Build Preconditioner 
        !$omp parallel do private(m,ibl,pl,i,j,M_block) schedule(guided)
        do m = 1, ns
          ibl = e3d%ib(m,1:27)
          pl  = e3d%pix(ibl(1:27))
		  
		  ! Fill the FULL 3x3 diagonal block for node 'm'
          do i = 1, 3
            do j = 1, 3
              M_block(i,j) = e3d%dk(pl(27),1,i,1,j) + e3d%dk(pl(7),2,i,2,j) +  &
                             e3d%dk(pl(6),3,i,3,j)  + e3d%dk(pl(5),4,i,4,j) +  &
                             e3d%dk(pl(25),5,i,5,j) + e3d%dk(pl(15),6,i,6,j) + &
                             e3d%dk(pl(14),7,i,7,j) + e3d%dk(pl(13),8,i,8,j)
            end do
          end do
		  
		  ! --- Invert block for each node ---
		  ! Stores the inverse to apply it rapidly during the CG iterations.
          call invert3x3(M_block, M_block_inv(m,:,:))
        end do
        !$omp end parallel do
      end if
  
      ! ======================================================================
      ! ---- 8. Preconditioned CG Algorithm Initialization -------------------
      ! Apply the inverse preconditioner to the initial residual (gradient).
      ! Mathematically: $z = M^{-1} * g$
      ! ======================================================================
      !$omp parallel do collapse(2) schedule(guided)
      do i = 1, ndof 
	    do m = 1, ns 
        z(m,i) = M_block_inv(m,i,1) * gb(m,1) + &
		         M_block_inv(m,i,2) * gb(m,2) + &
				 M_block_inv(m,i,3) * gb(m,3)
        end do
	  end do
      !$omp end parallel do
  
  
      ! ======================================================================
      ! Compute the preconditioned inner product: $gg = z^T * g$
      ! ======================================================================
      gg = 0.d0 
      !$omp parallel do collapse(2) schedule(guided) reduction(+:gg)	  
      do j = 1, ndof 
	    do i = 1, ns
          gg = gg + z(i,j)*gb(i,j)
        end do
	  end do
      !$omp end parallel do    

      ! ======================================================================
      ! Set the initial search direction $h$ equal to the preconditioned 
      ! residual $z$ on the first load step.
      ! ======================================================================
      if (kkk == 1) then
        !$omp parallel do collapse(2) schedule(guided)
        do i = 1, ndof 
		  do m = 1, ns
            h(m,i) = z(m,i)
          end do
        end do
        !$omp end parallel do
      end if
	  
      ! Initialize divergence tracking variables
      rel_res_prev = sqrt(gg/gg0)
      rel_res_min = rel_res_prev
      reset_count = 0	  


      ! ======================================================================
      ! ---- 9. Main Preconditioned Conjugate Gradient Loop ------------------
      ! ======================================================================
      Lstep = 0
      do ijk = 1, ldemb
        Lstep = Lstep + 1
		
		! Zero out the Matrix-Vector product array (Ah) before accumulation		
        Ah = 0.d0

        ! ====================================================================
        ! ---- 9a. Sparse Matrix-Vector Multiply (SpMV): Ah = A * h ----------
        ! Calculates the structural response to the current search direction.
        ! Processed in cache-friendly blocks to optimize memory bandwidth.
        ! ====================================================================
        do L1 = 1, ns_local, block_size
          L2 = min(L1 + block_size - 1, ns_local)

          ! Matrix multiply: Ah = A * h
          !$omp parallel do private(m, ibl, pl, ahsum1, ahsum2, ahsum3) schedule(guided)
          do m = L1, L2
            ! Caching node 'm' neighborhood connectivity and phase identifiers
            ibl = e3d%ib(m,1:27)      ! All 27 neighbors
            pl  = e3d%pix(ibl(1:27))  ! Phase index for all neighbors

            ! DOF 1 (X, j=1)
            ahsum1 = 0.d0
            ! -- n=1 --
            ahsum1 = ahsum1 &
          + h(ibl(1),1 )*( e3d%dk(pl(27),1,1,4,1) + e3d%dk(pl(7),2,1,3,1) + e3d%dk(pl(25),5,1,8,1) + e3d%dk(pl(15),6,1,7,1) ) &
          + h(ibl(2),1)*( e3d%dk(pl(27),1,1,3,1) + e3d%dk(pl(25),5,1,7,1) ) &
          + h(ibl(3),1)*( e3d%dk(pl(27),1,1,2,1) + e3d%dk(pl(5),4,1,3,1) + e3d%dk(pl(13),8,1,7,1) + e3d%dk(pl(25),5,1,6,1) ) &
          + h(ibl(4),1)*( e3d%dk(pl(5),4,1,2,1) + e3d%dk(pl(13),8,1,6,1) ) &
          + h(ibl(5),1)*( e3d%dk(pl(6),3,1,2,1) + e3d%dk(pl(5),4,1,1,1) + e3d%dk(pl(14),7,1,6,1) + e3d%dk(pl(13),8,1,5,1) ) &
          + h(ibl(6),1)*( e3d%dk(pl(6),3,1,1,1) + e3d%dk(pl(14),7,1,5,1) ) &
          + h(ibl(7),1)*( e3d%dk(pl(6),3,1,4,1) + e3d%dk(pl(7),2,1,1,1) + e3d%dk(pl(14),7,1,8,1) + e3d%dk(pl(15),6,1,5,1) ) &
          + h(ibl(8),1)*( e3d%dk(pl(7),2,1,4,1) + e3d%dk(pl(15),6,1,8,1) ) &
          + h(ibl(9),1)*( e3d%dk(pl(25),5,1,4,1) + e3d%dk(pl(15),6,1,3,1) ) &
          + h(ibl(10),1)*( e3d%dk(pl(25),5,1,3,1) ) &
          + h(ibl(11),1)*( e3d%dk(pl(13),8,1,3,1) + e3d%dk(pl(25),5,1,2,1) ) &
          + h(ibl(12),1)*( e3d%dk(pl(13),8,1,2,1) ) &
          + h(ibl(13),1)*( e3d%dk(pl(13),8,1,1,1) + e3d%dk(pl(14),7,1,2,1) ) &
          + h(ibl(14),1)*( e3d%dk(pl(14),7,1,1,1) ) &
          + h(ibl(15),1)*( e3d%dk(pl(14),7,1,4,1) + e3d%dk(pl(15),6,1,1,1) ) &
          + h(ibl(16),1)*( e3d%dk(pl(15),6,1,4,1) ) &
          + h(ibl(17),1)*( e3d%dk(pl(27),1,1,8,1) + e3d%dk(pl(7),2,1,7,1) ) &
          + h(ibl(18),1)*( e3d%dk(pl(27),1,1,7,1) ) &
          + h(ibl(19),1)*( e3d%dk(pl(27),1,1,6,1) + e3d%dk(pl(5),4,1,7,1) ) &
          + h(ibl(20),1)*( e3d%dk(pl(5),4,1,6,1) ) &
          + h(ibl(21),1)*( e3d%dk(pl(5),4,1,5,1) + e3d%dk(pl(6),3,1,6,1) ) &
          + h(ibl(22),1)*( e3d%dk(pl(6),3,1,5,1) ) &
          + h(ibl(23),1)*( e3d%dk(pl(6),3,1,8,1) + e3d%dk(pl(7),2,1,5,1) ) &
          + h(ibl(24),1)*( e3d%dk(pl(7),2,1,8,1) ) &
          + h(ibl(25),1)*( e3d%dk(pl(14),7,1,3,1) + e3d%dk(pl(13),8,1,4,1) + e3d%dk(pl(15),6,1,2,1) + e3d%dk(pl(25),5,1,1,1) ) &
          + h(ibl(26),1)*( e3d%dk(pl(6),3,1,7,1) + e3d%dk(pl(5),4,1,8,1) + e3d%dk(pl(27),1,1,5,1) + e3d%dk(pl(7),2,1,6,1) ) &
          + h(ibl(27),1)*( e3d%dk(pl(27),1,1,1,1) + e3d%dk(pl(7),2,1,2,1) + e3d%dk(pl(6),3,1,3,1) + e3d%dk(pl(5),4,1,4,1) &
                          + e3d%dk(pl(25),5,1,5,1) + e3d%dk(pl(15),6,1,6,1) + e3d%dk(pl(14),7,1,7,1) + e3d%dk(pl(13),8,1,8,1) )
            ! -- n=2 --
            ahsum1 = ahsum1 &
          + h(ibl(1),2 )*( e3d%dk(pl(27),1,1,4,2) + e3d%dk(pl(7),2,1,3,2) + e3d%dk(pl(25),5,1,8,2) + e3d%dk(pl(15),6,1,7,2) ) &
          + h(ibl(2),2)*( e3d%dk(pl(27),1,1,3,2) + e3d%dk(pl(25),5,1,7,2) ) &
          + h(ibl(3),2)*( e3d%dk(pl(27),1,1,2,2) + e3d%dk(pl(5),4,1,3,2) + e3d%dk(pl(13),8,1,7,2) + e3d%dk(pl(25),5,1,6,2) ) &
          + h(ibl(4),2)*( e3d%dk(pl(5),4,1,2,2) + e3d%dk(pl(13),8,1,6,2) ) &
          + h(ibl(5),2)*( e3d%dk(pl(6),3,1,2,2) + e3d%dk(pl(5),4,1,1,2) + e3d%dk(pl(14),7,1,6,2) + e3d%dk(pl(13),8,1,5,2) ) &
          + h(ibl(6),2)*( e3d%dk(pl(6),3,1,1,2) + e3d%dk(pl(14),7,1,5,2) ) &
          + h(ibl(7),2)*( e3d%dk(pl(6),3,1,4,2) + e3d%dk(pl(7),2,1,1,2) + e3d%dk(pl(14),7,1,8,2) + e3d%dk(pl(15),6,1,5,2) ) &
          + h(ibl(8),2)*( e3d%dk(pl(7),2,1,4,2) + e3d%dk(pl(15),6,1,8,2) ) &
          + h(ibl(9),2)*( e3d%dk(pl(25),5,1,4,2) + e3d%dk(pl(15),6,1,3,2) ) &
          + h(ibl(10),2)*( e3d%dk(pl(25),5,1,3,2) ) &
          + h(ibl(11),2)*( e3d%dk(pl(13),8,1,3,2) + e3d%dk(pl(25),5,1,2,2) ) &
          + h(ibl(12),2)*( e3d%dk(pl(13),8,1,2,2) ) &
          + h(ibl(13),2)*( e3d%dk(pl(13),8,1,1,2) + e3d%dk(pl(14),7,1,2,2) ) &
          + h(ibl(14),2)*( e3d%dk(pl(14),7,1,1,2) ) &
          + h(ibl(15),2)*( e3d%dk(pl(14),7,1,4,2) + e3d%dk(pl(15),6,1,1,2) ) &
          + h(ibl(16),2)*( e3d%dk(pl(15),6,1,4,2) ) &
          + h(ibl(17),2)*( e3d%dk(pl(27),1,1,8,2) + e3d%dk(pl(7),2,1,7,2) ) &
          + h(ibl(18),2)*( e3d%dk(pl(27),1,1,7,2) ) &
          + h(ibl(19),2)*( e3d%dk(pl(27),1,1,6,2) + e3d%dk(pl(5),4,1,7,2) ) &
          + h(ibl(20),2)*( e3d%dk(pl(5),4,1,6,2) ) &
          + h(ibl(21),2)*( e3d%dk(pl(5),4,1,5,2) + e3d%dk(pl(6),3,1,6,2) ) &
          + h(ibl(22),2)*( e3d%dk(pl(6),3,1,5,2) ) &
          + h(ibl(23),2)*( e3d%dk(pl(6),3,1,8,2) + e3d%dk(pl(7),2,1,5,2) ) &
          + h(ibl(24),2)*( e3d%dk(pl(7),2,1,8,2) ) &
          + h(ibl(25),2)*( e3d%dk(pl(14),7,1,3,2) + e3d%dk(pl(13),8,1,4,2) + e3d%dk(pl(15),6,1,2,2) + e3d%dk(pl(25),5,1,1,2) ) &
          + h(ibl(26),2)*( e3d%dk(pl(6),3,1,7,2) + e3d%dk(pl(5),4,1,8,2) + e3d%dk(pl(27),1,1,5,2) + e3d%dk(pl(7),2,1,6,2) ) &
          + h(ibl(27),2)*( e3d%dk(pl(27),1,1,1,2) + e3d%dk(pl(7),2,1,2,2) + e3d%dk(pl(6),3,1,3,2) + e3d%dk(pl(5),4,1,4,2) &
                          + e3d%dk(pl(25),5,1,5,2) + e3d%dk(pl(15),6,1,6,2) + e3d%dk(pl(14),7,1,7,2) + e3d%dk(pl(13),8,1,8,2) )
            ! -- n=3 --
            ahsum1 = ahsum1 &
          + h(ibl(1),3 )*( e3d%dk(pl(27),1,1,4,3) + e3d%dk(pl(7),2,1,3,3) + e3d%dk(pl(25),5,1,8,3) + e3d%dk(pl(15),6,1,7,3) ) &
          + h(ibl(2),3)*( e3d%dk(pl(27),1,1,3,3) + e3d%dk(pl(25),5,1,7,3) ) &
          + h(ibl(3),3)*( e3d%dk(pl(27),1,1,2,3) + e3d%dk(pl(5),4,1,3,3) + e3d%dk(pl(13),8,1,7,3) + e3d%dk(pl(25),5,1,6,3) ) &
          + h(ibl(4),3)*( e3d%dk(pl(5),4,1,2,3) + e3d%dk(pl(13),8,1,6,3) ) &
          + h(ibl(5),3)*( e3d%dk(pl(6),3,1,2,3) + e3d%dk(pl(5),4,1,1,3) + e3d%dk(pl(14),7,1,6,3) + e3d%dk(pl(13),8,1,5,3) ) &
          + h(ibl(6),3)*( e3d%dk(pl(6),3,1,1,3) + e3d%dk(pl(14),7,1,5,3) ) &
          + h(ibl(7),3)*( e3d%dk(pl(6),3,1,4,3) + e3d%dk(pl(7),2,1,1,3) + e3d%dk(pl(14),7,1,8,3) + e3d%dk(pl(15),6,1,5,3) ) &
          + h(ibl(8),3)*( e3d%dk(pl(7),2,1,4,3) + e3d%dk(pl(15),6,1,8,3) ) &
          + h(ibl(9),3)*( e3d%dk(pl(25),5,1,4,3) + e3d%dk(pl(15),6,1,3,3) ) &
          + h(ibl(10),3)*( e3d%dk(pl(25),5,1,3,3) ) &
          + h(ibl(11),3)*( e3d%dk(pl(13),8,1,3,3) + e3d%dk(pl(25),5,1,2,3) ) &
          + h(ibl(12),3)*( e3d%dk(pl(13),8,1,2,3) ) &
          + h(ibl(13),3)*( e3d%dk(pl(13),8,1,1,3) + e3d%dk(pl(14),7,1,2,3) ) &
          + h(ibl(14),3)*( e3d%dk(pl(14),7,1,1,3) ) &
          + h(ibl(15),3)*( e3d%dk(pl(14),7,1,4,3) + e3d%dk(pl(15),6,1,1,3) ) &
          + h(ibl(16),3)*( e3d%dk(pl(15),6,1,4,3) ) &
          + h(ibl(17),3)*( e3d%dk(pl(27),1,1,8,3) + e3d%dk(pl(7),2,1,7,3) ) &
          + h(ibl(18),3)*( e3d%dk(pl(27),1,1,7,3) ) &
          + h(ibl(19),3)*( e3d%dk(pl(27),1,1,6,3) + e3d%dk(pl(5),4,1,7,3) ) &
          + h(ibl(20),3)*( e3d%dk(pl(5),4,1,6,3) ) &
          + h(ibl(21),3)*( e3d%dk(pl(5),4,1,5,3) + e3d%dk(pl(6),3,1,6,3) ) &
          + h(ibl(22),3)*( e3d%dk(pl(6),3,1,5,3) ) &
          + h(ibl(23),3)*( e3d%dk(pl(6),3,1,8,3) + e3d%dk(pl(7),2,1,5,3) ) &
          + h(ibl(24),3)*( e3d%dk(pl(7),2,1,8,3) ) &
          + h(ibl(25),3)*( e3d%dk(pl(14),7,1,3,3) + e3d%dk(pl(13),8,1,4,3) + e3d%dk(pl(15),6,1,2,3) + e3d%dk(pl(25),5,1,1,3) ) &
          + h(ibl(26),3)*( e3d%dk(pl(6),3,1,7,3) + e3d%dk(pl(5),4,1,8,3) + e3d%dk(pl(27),1,1,5,3) + e3d%dk(pl(7),2,1,6,3) ) &
          + h(ibl(27),3)*( e3d%dk(pl(27),1,1,1,3) + e3d%dk(pl(7),2,1,2,3) + e3d%dk(pl(6),3,1,3,3) + e3d%dk(pl(5),4,1,4,3) &
                          + e3d%dk(pl(25),5,1,5,3) + e3d%dk(pl(15),6,1,6,3) + e3d%dk(pl(14),7,1,7,3) + e3d%dk(pl(13),8,1,8,3) )

            Ah(m,1) = ahsum1

            ! DOF 2 (Y, j=2)
            ahsum2 = 0.d0
            ! -- n=1 --
            ahsum2 = ahsum2 &
          + h(ibl(1),1 )*( e3d%dk(pl(27),1,2,4,1) + e3d%dk(pl(7),2,2,3,1) + e3d%dk(pl(25),5,2,8,1) + e3d%dk(pl(15),6,2,7,1) ) &
          + h(ibl(2),1)*( e3d%dk(pl(27),1,2,3,1) + e3d%dk(pl(25),5,2,7,1) ) &
          + h(ibl(3),1)*( e3d%dk(pl(27),1,2,2,1) + e3d%dk(pl(5),4,2,3,1) + e3d%dk(pl(13),8,2,7,1) + e3d%dk(pl(25),5,2,6,1) ) &
          + h(ibl(4),1)*( e3d%dk(pl(5),4,2,2,1) + e3d%dk(pl(13),8,2,6,1) ) &
          + h(ibl(5),1)*( e3d%dk(pl(6),3,2,2,1) + e3d%dk(pl(5),4,2,1,1) + e3d%dk(pl(14),7,2,6,1) + e3d%dk(pl(13),8,2,5,1) ) &
          + h(ibl(6),1)*( e3d%dk(pl(6),3,2,1,1) + e3d%dk(pl(14),7,2,5,1) ) &
          + h(ibl(7),1)*( e3d%dk(pl(6),3,2,4,1) + e3d%dk(pl(7),2,2,1,1) + e3d%dk(pl(14),7,2,8,1) + e3d%dk(pl(15),6,2,5,1) ) &
          + h(ibl(8),1)*( e3d%dk(pl(7),2,2,4,1) + e3d%dk(pl(15),6,2,8,1) ) &
          + h(ibl(9),1)*( e3d%dk(pl(25),5,2,4,1) + e3d%dk(pl(15),6,2,3,1) ) &
          + h(ibl(10),1)*( e3d%dk(pl(25),5,2,3,1) ) &
          + h(ibl(11),1)*( e3d%dk(pl(13),8,2,3,1) + e3d%dk(pl(25),5,2,2,1) ) &
          + h(ibl(12),1)*( e3d%dk(pl(13),8,2,2,1) ) &
          + h(ibl(13),1)*( e3d%dk(pl(13),8,2,1,1) + e3d%dk(pl(14),7,2,2,1) ) &
          + h(ibl(14),1)*( e3d%dk(pl(14),7,2,1,1) ) &
          + h(ibl(15),1)*( e3d%dk(pl(14),7,2,4,1) + e3d%dk(pl(15),6,2,1,1) ) &
          + h(ibl(16),1)*( e3d%dk(pl(15),6,2,4,1) ) &
          + h(ibl(17),1)*( e3d%dk(pl(27),1,2,8,1) + e3d%dk(pl(7),2,2,7,1) ) &
          + h(ibl(18),1)*( e3d%dk(pl(27),1,2,7,1) ) &
          + h(ibl(19),1)*( e3d%dk(pl(27),1,2,6,1) + e3d%dk(pl(5),4,2,7,1) ) &
          + h(ibl(20),1)*( e3d%dk(pl(5),4,2,6,1) ) &
          + h(ibl(21),1)*( e3d%dk(pl(5),4,2,5,1) + e3d%dk(pl(6),3,2,6,1) ) &
          + h(ibl(22),1)*( e3d%dk(pl(6),3,2,5,1) ) &
          + h(ibl(23),1)*( e3d%dk(pl(6),3,2,8,1) + e3d%dk(pl(7),2,2,5,1) ) &
          + h(ibl(24),1)*( e3d%dk(pl(7),2,2,8,1) ) &
          + h(ibl(25),1)*( e3d%dk(pl(14),7,2,3,1) + e3d%dk(pl(13),8,2,4,1) + e3d%dk(pl(15),6,2,2,1) + e3d%dk(pl(25),5,2,1,1) ) &
          + h(ibl(26),1)*( e3d%dk(pl(6),3,2,7,1) + e3d%dk(pl(5),4,2,8,1) + e3d%dk(pl(27),1,2,5,1) + e3d%dk(pl(7),2,2,6,1) ) &
          + h(ibl(27),1)*( e3d%dk(pl(27),1,2,1,1) + e3d%dk(pl(7),2,2,2,1) + e3d%dk(pl(6),3,2,3,1) + e3d%dk(pl(5),4,2,4,1) &
                          + e3d%dk(pl(25),5,2,5,1) + e3d%dk(pl(15),6,2,6,1) + e3d%dk(pl(14),7,2,7,1) + e3d%dk(pl(13),8,2,8,1) )
            ! -- n=2 --
            ahsum2 = ahsum2 &
          + h(ibl(1),2 )*( e3d%dk(pl(27),1,2,4,2) + e3d%dk(pl(7),2,2,3,2) + e3d%dk(pl(25),5,2,8,2) + e3d%dk(pl(15),6,2,7,2) ) &
          + h(ibl(2),2)*( e3d%dk(pl(27),1,2,3,2) + e3d%dk(pl(25),5,2,7,2) ) &
          + h(ibl(3),2)*( e3d%dk(pl(27),1,2,2,2) + e3d%dk(pl(5),4,2,3,2) + e3d%dk(pl(13),8,2,7,2) + e3d%dk(pl(25),5,2,6,2) ) &
          + h(ibl(4),2)*( e3d%dk(pl(5),4,2,2,2) + e3d%dk(pl(13),8,2,6,2) ) &
          + h(ibl(5),2)*( e3d%dk(pl(6),3,2,2,2) + e3d%dk(pl(5),4,2,1,2) + e3d%dk(pl(14),7,2,6,2) + e3d%dk(pl(13),8,2,5,2) ) &
          + h(ibl(6),2)*( e3d%dk(pl(6),3,2,1,2) + e3d%dk(pl(14),7,2,5,2) ) &
          + h(ibl(7),2)*( e3d%dk(pl(6),3,2,4,2) + e3d%dk(pl(7),2,2,1,2) + e3d%dk(pl(14),7,2,8,2) + e3d%dk(pl(15),6,2,5,2) ) &
          + h(ibl(8),2)*( e3d%dk(pl(7),2,2,4,2) + e3d%dk(pl(15),6,2,8,2) ) &
          + h(ibl(9),2)*( e3d%dk(pl(25),5,2,4,2) + e3d%dk(pl(15),6,2,3,2) ) &
          + h(ibl(10),2)*( e3d%dk(pl(25),5,2,3,2) ) &
          + h(ibl(11),2)*( e3d%dk(pl(13),8,2,3,2) + e3d%dk(pl(25),5,2,2,2) ) &
          + h(ibl(12),2)*( e3d%dk(pl(13),8,2,2,2) ) &
          + h(ibl(13),2)*( e3d%dk(pl(13),8,2,1,2) + e3d%dk(pl(14),7,2,2,2) ) &
          + h(ibl(14),2)*( e3d%dk(pl(14),7,2,1,2) ) &
          + h(ibl(15),2)*( e3d%dk(pl(14),7,2,4,2) + e3d%dk(pl(15),6,2,1,2) ) &
          + h(ibl(16),2)*( e3d%dk(pl(15),6,2,4,2) ) &
          + h(ibl(17),2)*( e3d%dk(pl(27),1,2,8,2) + e3d%dk(pl(7),2,2,7,2) ) &
          + h(ibl(18),2)*( e3d%dk(pl(27),1,2,7,2) ) &
          + h(ibl(19),2)*( e3d%dk(pl(27),1,2,6,2) + e3d%dk(pl(5),4,2,7,2) ) &
          + h(ibl(20),2)*( e3d%dk(pl(5),4,2,6,2) ) &
          + h(ibl(21),2)*( e3d%dk(pl(5),4,2,5,2) + e3d%dk(pl(6),3,2,6,2) ) &
          + h(ibl(22),2)*( e3d%dk(pl(6),3,2,5,2) ) &
          + h(ibl(23),2)*( e3d%dk(pl(6),3,2,8,2) + e3d%dk(pl(7),2,2,5,2) ) &
          + h(ibl(24),2)*( e3d%dk(pl(7),2,2,8,2) ) &
          + h(ibl(25),2)*( e3d%dk(pl(14),7,2,3,2) + e3d%dk(pl(13),8,2,4,2) + e3d%dk(pl(15),6,2,2,2) + e3d%dk(pl(25),5,2,1,2) ) &
          + h(ibl(26),2)*( e3d%dk(pl(6),3,2,7,2) + e3d%dk(pl(5),4,2,8,2) + e3d%dk(pl(27),1,2,5,2) + e3d%dk(pl(7),2,2,6,2) ) &
          + h(ibl(27),2)*( e3d%dk(pl(27),1,2,1,2) + e3d%dk(pl(7),2,2,2,2) + e3d%dk(pl(6),3,2,3,2) + e3d%dk(pl(5),4,2,4,2) &
                          + e3d%dk(pl(25),5,2,5,2) + e3d%dk(pl(15),6,2,6,2) + e3d%dk(pl(14),7,2,7,2) + e3d%dk(pl(13),8,2,8,2) )
            ! -- n=3 --
            ahsum2 = ahsum2 &
          + h(ibl(1),3 )*( e3d%dk(pl(27),1,2,4,3) + e3d%dk(pl(7),2,2,3,3) + e3d%dk(pl(25),5,2,8,3) + e3d%dk(pl(15),6,2,7,3) ) &
          + h(ibl(2),3)*( e3d%dk(pl(27),1,2,3,3) + e3d%dk(pl(25),5,2,7,3) ) &
          + h(ibl(3),3)*( e3d%dk(pl(27),1,2,2,3) + e3d%dk(pl(5),4,2,3,3) + e3d%dk(pl(13),8,2,7,3) + e3d%dk(pl(25),5,2,6,3) ) &
          + h(ibl(4),3)*( e3d%dk(pl(5),4,2,2,3) + e3d%dk(pl(13),8,2,6,3) ) &
          + h(ibl(5),3)*( e3d%dk(pl(6),3,2,2,3) + e3d%dk(pl(5),4,2,1,3) + e3d%dk(pl(14),7,2,6,3) + e3d%dk(pl(13),8,2,5,3) ) &
          + h(ibl(6),3)*( e3d%dk(pl(6),3,2,1,3) + e3d%dk(pl(14),7,2,5,3) ) &
          + h(ibl(7),3)*( e3d%dk(pl(6),3,2,4,3) + e3d%dk(pl(7),2,2,1,3) + e3d%dk(pl(14),7,2,8,3) + e3d%dk(pl(15),6,2,5,3) ) &
          + h(ibl(8),3)*( e3d%dk(pl(7),2,2,4,3) + e3d%dk(pl(15),6,2,8,3) ) &
          + h(ibl(9),3)*( e3d%dk(pl(25),5,2,4,3) + e3d%dk(pl(15),6,2,3,3) ) &
          + h(ibl(10),3)*( e3d%dk(pl(25),5,2,3,3) ) &
          + h(ibl(11),3)*( e3d%dk(pl(13),8,2,3,3) + e3d%dk(pl(25),5,2,2,3) ) &
          + h(ibl(12),3)*( e3d%dk(pl(13),8,2,2,3) ) &
          + h(ibl(13),3)*( e3d%dk(pl(13),8,2,1,3) + e3d%dk(pl(14),7,2,2,3) ) &
          + h(ibl(14),3)*( e3d%dk(pl(14),7,2,1,3) ) &
          + h(ibl(15),3)*( e3d%dk(pl(14),7,2,4,3) + e3d%dk(pl(15),6,2,1,3) ) &
          + h(ibl(16),3)*( e3d%dk(pl(15),6,2,4,3) ) &
          + h(ibl(17),3)*( e3d%dk(pl(27),1,2,8,3) + e3d%dk(pl(7),2,2,7,3) ) &
          + h(ibl(18),3)*( e3d%dk(pl(27),1,2,7,3) ) &
          + h(ibl(19),3)*( e3d%dk(pl(27),1,2,6,3) + e3d%dk(pl(5),4,2,7,3) ) &
          + h(ibl(20),3)*( e3d%dk(pl(5),4,2,6,3) ) &
          + h(ibl(21),3)*( e3d%dk(pl(5),4,2,5,3) + e3d%dk(pl(6),3,2,6,3) ) &
          + h(ibl(22),3)*( e3d%dk(pl(6),3,2,5,3) ) &
          + h(ibl(23),3)*( e3d%dk(pl(6),3,2,8,3) + e3d%dk(pl(7),2,2,5,3) ) &
          + h(ibl(24),3)*( e3d%dk(pl(7),2,2,8,3) ) &
          + h(ibl(25),3)*( e3d%dk(pl(14),7,2,3,3) + e3d%dk(pl(13),8,2,4,3) + e3d%dk(pl(15),6,2,2,3) + e3d%dk(pl(25),5,2,1,3) ) &
          + h(ibl(26),3)*( e3d%dk(pl(6),3,2,7,3) + e3d%dk(pl(5),4,2,8,3) + e3d%dk(pl(27),1,2,5,3) + e3d%dk(pl(7),2,2,6,3) ) &
          + h(ibl(27),3)*( e3d%dk(pl(27),1,2,1,3) + e3d%dk(pl(7),2,2,2,3) + e3d%dk(pl(6),3,2,3,3) + e3d%dk(pl(5),4,2,4,3) &
                          + e3d%dk(pl(25),5,2,5,3) + e3d%dk(pl(15),6,2,6,3) + e3d%dk(pl(14),7,2,7,3) + e3d%dk(pl(13),8,2,8,3) )

            Ah(m,2) = ahsum2

            ! DOF 3 (Z, j=3)
            ahsum3 = 0.d0
            ! -- n=1 --
            ahsum3 = ahsum3 &
          + h(ibl(1),1 )*( e3d%dk(pl(27),1,3,4,1) + e3d%dk(pl(7),2,3,3,1) + e3d%dk(pl(25),5,3,8,1) + e3d%dk(pl(15),6,3,7,1) ) &
          + h(ibl(2),1)*( e3d%dk(pl(27),1,3,3,1) + e3d%dk(pl(25),5,3,7,1) ) &
          + h(ibl(3),1)*( e3d%dk(pl(27),1,3,2,1) + e3d%dk(pl(5),4,3,3,1) + e3d%dk(pl(13),8,3,7,1) + e3d%dk(pl(25),5,3,6,1) ) &
          + h(ibl(4),1)*( e3d%dk(pl(5),4,3,2,1) + e3d%dk(pl(13),8,3,6,1) ) &
          + h(ibl(5),1)*( e3d%dk(pl(6),3,3,2,1) + e3d%dk(pl(5),4,3,1,1) + e3d%dk(pl(14),7,3,6,1) + e3d%dk(pl(13),8,3,5,1) ) &
          + h(ibl(6),1)*( e3d%dk(pl(6),3,3,1,1) + e3d%dk(pl(14),7,3,5,1) ) &
          + h(ibl(7),1)*( e3d%dk(pl(6),3,3,4,1) + e3d%dk(pl(7),2,3,1,1) + e3d%dk(pl(14),7,3,8,1) + e3d%dk(pl(15),6,3,5,1) ) &
          + h(ibl(8),1)*( e3d%dk(pl(7),2,3,4,1) + e3d%dk(pl(15),6,3,8,1) ) &
          + h(ibl(9),1)*( e3d%dk(pl(25),5,3,4,1) + e3d%dk(pl(15),6,3,3,1) ) &
          + h(ibl(10),1)*( e3d%dk(pl(25),5,3,3,1) ) &
          + h(ibl(11),1)*( e3d%dk(pl(13),8,3,3,1) + e3d%dk(pl(25),5,3,2,1) ) &
          + h(ibl(12),1)*( e3d%dk(pl(13),8,3,2,1) ) &
          + h(ibl(13),1)*( e3d%dk(pl(13),8,3,1,1) + e3d%dk(pl(14),7,3,2,1) ) &
          + h(ibl(14),1)*( e3d%dk(pl(14),7,3,1,1) ) &
          + h(ibl(15),1)*( e3d%dk(pl(14),7,3,4,1) + e3d%dk(pl(15),6,3,1,1) ) &
          + h(ibl(16),1)*( e3d%dk(pl(15),6,3,4,1) ) &
          + h(ibl(17),1)*( e3d%dk(pl(27),1,3,8,1) + e3d%dk(pl(7),2,3,7,1) ) &
          + h(ibl(18),1)*( e3d%dk(pl(27),1,3,7,1) ) &
          + h(ibl(19),1)*( e3d%dk(pl(27),1,3,6,1) + e3d%dk(pl(5),4,3,7,1) ) &
          + h(ibl(20),1)*( e3d%dk(pl(5),4,3,6,1) ) &
          + h(ibl(21),1)*( e3d%dk(pl(5),4,3,5,1) + e3d%dk(pl(6),3,3,6,1) ) &
          + h(ibl(22),1)*( e3d%dk(pl(6),3,3,5,1) ) &
          + h(ibl(23),1)*( e3d%dk(pl(6),3,3,8,1) + e3d%dk(pl(7),2,3,5,1) ) &
          + h(ibl(24),1)*( e3d%dk(pl(7),2,3,8,1) ) &
          + h(ibl(25),1)*( e3d%dk(pl(14),7,3,3,1) + e3d%dk(pl(13),8,3,4,1) + e3d%dk(pl(15),6,3,2,1) + e3d%dk(pl(25),5,3,1,1) ) &
          + h(ibl(26),1)*( e3d%dk(pl(6),3,3,7,1) + e3d%dk(pl(5),4,3,8,1) + e3d%dk(pl(27),1,3,5,1) + e3d%dk(pl(7),2,3,6,1) ) &
          + h(ibl(27),1)*( e3d%dk(pl(27),1,3,1,1) + e3d%dk(pl(7),2,3,2,1) + e3d%dk(pl(6),3,3,3,1) + e3d%dk(pl(5),4,3,4,1) &
                          + e3d%dk(pl(25),5,3,5,1) + e3d%dk(pl(15),6,3,6,1) + e3d%dk(pl(14),7,3,7,1) + e3d%dk(pl(13),8,3,8,1) )
            ! -- n=2 --
            ahsum3 = ahsum3 &
          + h(ibl(1),2 )*( e3d%dk(pl(27),1,3,4,2) + e3d%dk(pl(7),2,3,3,2) + e3d%dk(pl(25),5,3,8,2) + e3d%dk(pl(15),6,3,7,2) ) &
          + h(ibl(2),2)*( e3d%dk(pl(27),1,3,3,2) + e3d%dk(pl(25),5,3,7,2) ) &
          + h(ibl(3),2)*( e3d%dk(pl(27),1,3,2,2) + e3d%dk(pl(5),4,3,3,2) + e3d%dk(pl(13),8,3,7,2) + e3d%dk(pl(25),5,3,6,2) ) &
          + h(ibl(4),2)*( e3d%dk(pl(5),4,3,2,2) + e3d%dk(pl(13),8,3,6,2) ) &
          + h(ibl(5),2)*( e3d%dk(pl(6),3,3,2,2) + e3d%dk(pl(5),4,3,1,2) + e3d%dk(pl(14),7,3,6,2) + e3d%dk(pl(13),8,3,5,2) ) &
          + h(ibl(6),2)*( e3d%dk(pl(6),3,3,1,2) + e3d%dk(pl(14),7,3,5,2) ) &
          + h(ibl(7),2)*( e3d%dk(pl(6),3,3,4,2) + e3d%dk(pl(7),2,3,1,2) + e3d%dk(pl(14),7,3,8,2) + e3d%dk(pl(15),6,3,5,2) ) &
          + h(ibl(8),2)*( e3d%dk(pl(7),2,3,4,2) + e3d%dk(pl(15),6,3,8,2) ) &
          + h(ibl(9),2)*( e3d%dk(pl(25),5,3,4,2) + e3d%dk(pl(15),6,3,3,2) ) &
          + h(ibl(10),2)*( e3d%dk(pl(25),5,3,3,2) ) &
          + h(ibl(11),2)*( e3d%dk(pl(13),8,3,3,2) + e3d%dk(pl(25),5,3,2,2) ) &
          + h(ibl(12),2)*( e3d%dk(pl(13),8,3,2,2) ) &
          + h(ibl(13),2)*( e3d%dk(pl(13),8,3,1,2) + e3d%dk(pl(14),7,3,2,2) ) &
          + h(ibl(14),2)*( e3d%dk(pl(14),7,3,1,2) ) &
          + h(ibl(15),2)*( e3d%dk(pl(14),7,3,4,2) + e3d%dk(pl(15),6,3,1,2) ) &
          + h(ibl(16),2)*( e3d%dk(pl(15),6,3,4,2) ) &
          + h(ibl(17),2)*( e3d%dk(pl(27),1,3,8,2) + e3d%dk(pl(7),2,3,7,2) ) &
          + h(ibl(18),2)*( e3d%dk(pl(27),1,3,7,2) ) &
          + h(ibl(19),2)*( e3d%dk(pl(27),1,3,6,2) + e3d%dk(pl(5),4,3,7,2) ) &
          + h(ibl(20),2)*( e3d%dk(pl(5),4,3,6,2) ) &
          + h(ibl(21),2)*( e3d%dk(pl(5),4,3,5,2) + e3d%dk(pl(6),3,3,6,2) ) &
          + h(ibl(22),2)*( e3d%dk(pl(6),3,3,5,2) ) &
          + h(ibl(23),2)*( e3d%dk(pl(6),3,3,8,2) + e3d%dk(pl(7),2,3,5,2) ) &
          + h(ibl(24),2)*( e3d%dk(pl(7),2,3,8,2) ) &
          + h(ibl(25),2)*( e3d%dk(pl(14),7,3,3,2) + e3d%dk(pl(13),8,3,4,2) + e3d%dk(pl(15),6,3,2,2) + e3d%dk(pl(25),5,3,1,2) ) &
          + h(ibl(26),2)*( e3d%dk(pl(6),3,3,7,2) + e3d%dk(pl(5),4,3,8,2) + e3d%dk(pl(27),1,3,5,2) + e3d%dk(pl(7),2,3,6,2) ) &
          + h(ibl(27),2)*( e3d%dk(pl(27),1,3,1,2) + e3d%dk(pl(7),2,3,2,2) + e3d%dk(pl(6),3,3,3,2) + e3d%dk(pl(5),4,3,4,2) &
                          + e3d%dk(pl(25),5,3,5,2) + e3d%dk(pl(15),6,3,6,2) + e3d%dk(pl(14),7,3,7,2) + e3d%dk(pl(13),8,3,8,2) )
            ! -- n=3 --
            ahsum3 = ahsum3 &
          + h(ibl(1),3 )*( e3d%dk(pl(27),1,3,4,3) + e3d%dk(pl(7),2,3,3,3) + e3d%dk(pl(25),5,3,8,3) + e3d%dk(pl(15),6,3,7,3) ) &
          + h(ibl(2),3)*( e3d%dk(pl(27),1,3,3,3) + e3d%dk(pl(25),5,3,7,3) ) &
          + h(ibl(3),3)*( e3d%dk(pl(27),1,3,2,3) + e3d%dk(pl(5),4,3,3,3) + e3d%dk(pl(13),8,3,7,3) + e3d%dk(pl(25),5,3,6,3) ) &
          + h(ibl(4),3)*( e3d%dk(pl(5),4,3,2,3) + e3d%dk(pl(13),8,3,6,3) ) &
          + h(ibl(5),3)*( e3d%dk(pl(6),3,3,2,3) + e3d%dk(pl(5),4,3,1,3) + e3d%dk(pl(14),7,3,6,3) + e3d%dk(pl(13),8,3,5,3) ) &
          + h(ibl(6),3)*( e3d%dk(pl(6),3,3,1,3) + e3d%dk(pl(14),7,3,5,3) ) &
          + h(ibl(7),3)*( e3d%dk(pl(6),3,3,4,3) + e3d%dk(pl(7),2,3,1,3) + e3d%dk(pl(14),7,3,8,3) + e3d%dk(pl(15),6,3,5,3) ) &
          + h(ibl(8),3)*( e3d%dk(pl(7),2,3,4,3) + e3d%dk(pl(15),6,3,8,3) ) &
          + h(ibl(9),3)*( e3d%dk(pl(25),5,3,4,3) + e3d%dk(pl(15),6,3,3,3) ) &
          + h(ibl(10),3)*( e3d%dk(pl(25),5,3,3,3) ) &
          + h(ibl(11),3)*( e3d%dk(pl(13),8,3,3,3) + e3d%dk(pl(25),5,3,2,3) ) &
          + h(ibl(12),3)*( e3d%dk(pl(13),8,3,2,3) ) &
          + h(ibl(13),3)*( e3d%dk(pl(13),8,3,1,3) + e3d%dk(pl(14),7,3,2,3) ) &
          + h(ibl(14),3)*( e3d%dk(pl(14),7,3,1,3) ) &
          + h(ibl(15),3)*( e3d%dk(pl(14),7,3,4,3) + e3d%dk(pl(15),6,3,1,3) ) &
          + h(ibl(16),3)*( e3d%dk(pl(15),6,3,4,3) ) &
          + h(ibl(17),3)*( e3d%dk(pl(27),1,3,8,3) + e3d%dk(pl(7),2,3,7,3) ) &
          + h(ibl(18),3)*( e3d%dk(pl(27),1,3,7,3) ) &
          + h(ibl(19),3)*( e3d%dk(pl(27),1,3,6,3) + e3d%dk(pl(5),4,3,7,3) ) &
          + h(ibl(20),3)*( e3d%dk(pl(5),4,3,6,3) ) &
          + h(ibl(21),3)*( e3d%dk(pl(5),4,3,5,3) + e3d%dk(pl(6),3,3,6,3) ) &
          + h(ibl(22),3)*( e3d%dk(pl(6),3,3,5,3) ) &
          + h(ibl(23),3)*( e3d%dk(pl(6),3,3,8,3) + e3d%dk(pl(7),2,3,5,3) ) &
          + h(ibl(24),3)*( e3d%dk(pl(7),2,3,8,3) ) &
          + h(ibl(25),3)*( e3d%dk(pl(14),7,3,3,3) + e3d%dk(pl(13),8,3,4,3) + e3d%dk(pl(15),6,3,2,3) + e3d%dk(pl(25),5,3,1,3) ) &
          + h(ibl(26),3)*( e3d%dk(pl(6),3,3,7,3) + e3d%dk(pl(5),4,3,8,3) + e3d%dk(pl(27),1,3,5,3) + e3d%dk(pl(7),2,3,6,3) ) &
          + h(ibl(27),3)*( e3d%dk(pl(27),1,3,1,3) + e3d%dk(pl(7),2,3,2,3) + e3d%dk(pl(6),3,3,3,3) + e3d%dk(pl(5),4,3,4,3) &
                          + e3d%dk(pl(25),5,3,5,3) + e3d%dk(pl(15),6,3,6,3) + e3d%dk(pl(14),7,3,7,3) + e3d%dk(pl(13),8,3,8,3) )

            Ah(m,3) = ahsum3

          end do
          !$omp end parallel do

        end do ! blocks

        ! ====================================================================
        ! ---- 9b. Curvature and Step Size Calculation -----------------------
        ! Compute the curvature along the search direction: hAh = h^T * (A * h)
        ! ====================================================================
        hAh = 0.d0
        !$omp parallel do collapse(2) schedule(guided) reduction(+:hAh)		
        do i = 1, ndof
		  do m = 1, ns
            hAh = hAh + h(m,i) * Ah(m,i)
          end do
        end do
        !$omp end parallel do

        ! Protect against division by zero if the curvature vanishes
        if(dabs(hAh) < 1d-14) then
          print *,"WARNING: hAh is nearly zero, breaking CG step."
          exit
        end if
		
		
		! Compute Step Size (lambda, often denoted \alpha in literature)
        ! $\lambda = (z^T g) / (h^T A h)$
        lambda = gg / hAh

        ! ====================================================================
        ! ---- 9c. Solution and Residual Update ------------------------------
        ! Update the displacement field (u) and the negative gradient/residual (gb)
        ! using the calculated step size.
        ! ====================================================================
        !$omp parallel do collapse(2) schedule(guided)		
        do i = 1, ndof
		  do m = 1, ns
            u(m,i) = u(m,i) - lambda * h(m,i)
            gb(m,i) = gb(m,i) - lambda * Ah(m,i)
		  end do
        end do
        !$omp end parallel do

        ! ====================================================================
        ! ---- 9d. Preconditioner Application & New Residual -----------------
        ! Apply the Block Jacobi preconditioner to the newly updated residual.
        ! ====================================================================
        !$omp parallel do collapse(2) schedule(guided)		
        do i = 1, ndof
		  do m = 1, ns
            z(m,i) = M_block_inv(m,i,1)*gb(m,1) + &
		             M_block_inv(m,i,2)*gb(m,2) + &
			         M_block_inv(m,i,3)*gb(m,3)
          end do
        end do
        !$omp end parallel do
    
        ! Store previous inner product and compute the new one for the next step
        gglast = gg
        gg = 0.d0
        !$omp parallel do collapse(2) schedule(guided) reduction(+:gg)		
        do i = 1, ndof
		  do m = 1, ns
            gg = gg + z(m,i) * gb(m,i)
          end do
        end do
        !$omp end parallel do    

        ! Compute relative residual for convergence and divergence checks
		rel_res = sqrt(gg/gg0)
		
		! Store history of iterations
        history_step(Lstep) = Lstep
        history_res(Lstep)  = rel_res	

        ! ====================================================================
        ! ---- 9e. Divergence Control ----------------------------------------
        ! If the residual spikes (increases by an order of magnitude), the 
        ! conjugate directions have likely lost orthogonality due to numerical 
        ! noise. We wipe the history and fall back to Preconditioned Steepest Descent.
        ! ====================================================================
        if (rel_res > 10.0d0 * rel_res_prev .and. Lstep > 1) then
          print *, 'WARNING: Residual increased by more than 10x!'
          print *, '  Previous rel_res =', rel_res_prev
          print *, '  Current rel_res  =', rel_res
          print *, '  Resetting search direction to steepest descent.'
          
          ! Reset direction to preconditioned residual (steepest descent)
          !$omp parallel do collapse(2) schedule(guided)		
          do i = 1, ndof
            do m = 1, ns
              h(m,i) = z(m,i)
            end do
          end do
          !$omp end parallel do
          
          ! Reset gamma tracking
          rel_res_min = rel_res
          reset_count = reset_count + 1
          
          print *, '  Direction reset count:', reset_count
		  
        else
          ! ==================================================================
          ! ---- 9f. Conjugate Direction Update (Fletcher-Reeves) ------------
          ! Calculate orthogonalization factor (gamma, often denoted \beta).
          ! $\gamma = (z_{k+1}^T g_{k+1}) / (z_k^T g_k)$
          ! ==================================================================
          gamma = gg / gglast
          !$omp parallel do collapse(2) schedule(guided)		
          do i = 1, ndof
		    do m = 1, ns
		      h(m,i) = z(m,i) + gamma * h(m,i)
		    end do
          end do
          !$omp end parallel do     
		  
          ! Track minimum residual
          if (rel_res < rel_res_min) then
            rel_res_min = rel_res
          end if		  
        end if    

        ! Update previous residual for next iteration
        rel_res_prev = rel_res
        
        ! Exit early if converged (assuming 'tol' is defined globally/in module)		
        if (rel_res < tol) exit			
    
        print *, 'Number of conjugate steps_local iteration =', Lstep
        print *, 'rel. residual norm=', rel_res

      end do ! CG steps
	  
      ! ======================================================================
      ! ---- 10. File I/O & Global State Persistence -------------------------
      ! Write convergence history to disk and map local solver arrays back 
      ! into the global elasticity data structure (e3d).
      ! ======================================================================
      do i = 1, Lstep
         write(100, '(I8, 4X, E20.8)') history_step(i), history_res(i)
      end do

      ! Save back to struct fields	  
      !$omp parallel do collapse(2) schedule(guided)  
      do i = 1, ndof
	    do m = 1, ns
          e3d%u(m,i)  = u(m,i)
          e3d%gb(m,i) = gb(m,i)
          e3d%h(m,i)  = h(m,i)
        end do
      end do
      !$omp end parallel do	  


      ! ======================================================================
      ! ---- 11. Memory Cleanup ----------------------------------------------
      ! Prevent memory leaks upon returning to the main program.
      ! ======================================================================
      deallocate( u, gb, h, Ah )
      deallocate( M_block_inv, z )  
      deallocate( history_step, history_res )		  

      end subroutine dembx_OpenMP	  
	  
!-------------------------------------------------------------------------  
      pure subroutine invert3x3(a, ainv)
      ! ==============================================================================
      ! SUBROUTINE: invert3x3
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   A highly optimized, thread-safe utility to calculate the inverse of a 
      !   3x3 matrix. It is primarily used by the PCG solver to invert the local 
      !   Block Jacobi preconditioner matrices for each finite element node.
      !
      ! DUMMY ARGUMENTS:
      !   - a (intent(in)): The input 3x3 real matrix to be inverted.
      !   - ainv (intent(out)): The output 3x3 real inverted matrix.
      !
      ! METHODOLOGY:
      !   - Thread Safety: Declared as a 'pure' subroutine, guaranteeing no side 
      !     effects or global state modifications, making it completely safe to 
      !     call from within concurrent OpenMP loops.
      !   - Algorithm: Uses explicit cofactor expansion (Cramer's Rule / Adjugate 
      !     matrix method) to directly compute the inverse.
      !   - Optimization: Computes only the first column of cofactors initially 
      !     to find the determinant. If singular, it exits early to save operations. 
      !     It also computes the inverse determinant once (1 division) and uses 
      !     multiplication for the remaining 9 elements to maximize floating-point 
      !     throughput.
      !   - Stability: Includes a singularity check (abs(det) < 1d-14). If the 
      !     matrix is non-invertible, it gracefully returns a zero matrix.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none
	  
      real(dp), intent(in)  :: a(3,3)
      real(dp), intent(out) :: ainv(3,3)
	  
      real(dp) :: det, inv_det
	  ! Pre-calculate cofactors to save ops
      real(dp) :: c11, c12, c13
      real(dp) :: c21, c22, c23
      real(dp) :: c31, c32, c33

      ! ======================================================================
      ! ---- 1. Calculate Cofactors for Column 1 of Inverse ------------------
      ! These correspond to the Row 1 cofactors of the input matrix A.
      ! ======================================================================
      c11 =  a(2,2)*a(3,3) - a(2,3)*a(3,2)
      c21 = -a(2,1)*a(3,3) + a(2,3)*a(3,1)
      c31 =  a(2,1)*a(3,2) - a(2,2)*a(3,1)

      ! ======================================================================
      ! ---- 2. Compute Determinant ------------------------------------------
      ! Reuse c11, c21, and c31 to find the determinant via Laplace expansion.
      ! $Det(A) = a_{11}C_{11} + a_{12}C_{12} + a_{13}C_{13}$
      ! ======================================================================
      det = a(1,1)*c11 + a(1,2)*c21 + a(1,3)*c31
			

      ! ======================================================================
      ! ---- 3. Singularity Check --------------------------------------------
      ! If the determinant is near zero, return a null matrix to avoid 
      ! floating-point overflow during division.
      ! ======================================================================
      if (abs(det) < 1d-14) then
        ainv(:,:) = 0.d0
        return
      end if
	  
	  ! Pre-calculate the scalar reciprocal for the final scaling step
      inv_det = 1.0d0 / det	  
      
	  
      ! ======================================================================
      ! ---- 4. Calculate Remaining Cofactors --------------------------------
      ! These define the 2nd and 3rd columns of the inverse matrix.
      ! ====================================================================== 
      ! Column 2 of Inverse (Row 2 cofactors of A)
      c12 = -a(1,2)*a(3,3) + a(1,3)*a(3,2)
      c22 =  a(1,1)*a(3,3) - a(1,3)*a(3,1)
      c32 = -a(1,1)*a(3,2) + a(1,2)*a(3,1)

      ! Column 3 of Inverse (Row 3 cofactors of A)
      c13 =  a(1,2)*a(2,3) - a(1,3)*a(2,2)
      c23 = -a(1,1)*a(2,3) + a(1,3)*a(2,1)
      c33 =  a(1,1)*a(2,2) - a(1,2)*a(2,1)	  
		
		
      ! ======================================================================
      ! ---- 5. Final Assembly -----------------------------------------------
      ! Map cofactors to the output array and scale by the inverse determinant.
      ! The assignments follow column-major order for optimal memory writing.
      ! ======================================================================
      ainv(1,1) = c11 * inv_det
      ainv(2,1) = c21 * inv_det
      ainv(3,1) = c31 * inv_det

      ainv(1,2) = c12 * inv_det
      ainv(2,2) = c22 * inv_det
      ainv(3,2) = c32 * inv_det

      ainv(1,3) = c13 * inv_det
      ainv(2,3) = c23 * inv_det
      ainv(3,3) = c33 * inv_det

      end subroutine invert3x3		  
!-------------------------------------------------------------------------
      subroutine stress_OpenMP(e3d)
      ! ==============================================================================
      ! SUBROUTINE: stress_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Calculates the local and volume-averaged stress and strain tensors across
      !   the entire 3D domain. It post-processes the displacement field (u) by 
      !   computing gradients at voxel centers and applying Hooke's Law.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(inout)): Derived type containing global fields, including 
      !     solved displacements and the global average stress/strain outputs.
      !
      ! 'e3d' COMPONENTS MODIFIED/ACCESSED:
      !   - e3d%u (intent(in)): The solved nodal displacement field.
      !   - e3d%cmod (intent(in)): The 6x6 stiffness tensors for each material phase.
      !   - e3d%ib (intent(in)): Periodic connectivity mapping.
      !   - e3d%strxx... (intent(out)): Volume-averaged stress components.
      !   - e3d%sxx... (intent(out)): Volume-averaged strain components.
      !
      ! METHODOLOGY:
      !   - Finite Element: Uses 8-node trilinear hexahedral elements (voxels).
      !   - Strain-Displacement: Explicitly builds the B-matrix (es) using 
      !     pre-computed shape function derivatives.
      !   - Boundary Conditions: Re-applies the macroscopic strain jumps during 
      !     nodal displacement gathering to ensure consistency with periodic BCs.
      !   - Parallelization: Uses OpenMP over the 3D grid with `reduction` clauses 
      !     to safely compute global volume averages across multiple threads.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none

      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d

      ! Loop variables
      integer :: i, j, k, m,  mm, n8, n3, n
	  integer :: m_local, k_local, j_local, i_local, pixm

	  ! Geometry and Connectivity
      integer :: ibl(27)
      integer :: nxy
	  
      ! Element Strain-Displacement Matrix (6 components, 8 nodes, 3 DOFs)
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
	  real(dp) :: es(6,nnode_fe,ndof)
	  
	  ! Local Displacement and Results
	  real(dp) :: uu(nnode_fe,ndof)
      real(dp) :: str11, str22, str33, str13, str23, str12
      real(dp) :: s11, s22, s33, s13, s23, s12
	  
	  ! Pixel accumulators
      real(dp) :: str11_t, str22_t, str33_t, str13_t, str23_t, str12_t
      real(dp) :: s11_t, s22_t, s33_t, s13_t, s23_t, s12_t      
	  
	  ! Global Accumulators
      real(dp) :: strxx, stryy, strzz, strxz, stryz, strxy
      real(dp) :: sxx, syy, szz, sxz, syz, sxy
	  
      ! Macroscopic strains
	  real(dp) :: exx, eyy, ezz, exz, eyz, exy



      ! ======================================================================
      ! ---- 1. Initialization and B-Matrix Setup ----------------------------
      ! Unpack macroscopic strains and define shape function derivatives for 
      ! the center of the voxel.
      ! ======================================================================
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy
      nxy = nx * ny

      ! Pre-compute Shape Function Derivatives
      dndx = [-0.25d0, 0.25d0, 0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0, -0.25d0]
      dndy = [-0.25d0, -0.25d0, 0.25d0, 0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0]
      dndz = [-0.25d0, -0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0]

      ! Build Strain-Displacement Matrix (6 components x 8 nodes x 3 DOFs)
      ! Mapping: 1=xx, 2=yy, 3=zz, 4=xz, 5=yz, 6=xy
      es = 0.d0
      do n = 1, nnode_fe
        es(1,n,1) = dndx(n)
        es(2,n,2) = dndy(n)
        es(3,n,3) = dndz(n)
        es(4,n,1) = dndz(n);   es(4,n,3) = dndx(n); !  xz
        es(5,n,2) = dndz(n);   es(5,n,3) = dndy(n); !  yz
        es(6,n,1) = dndy(n);   es(6,n,2) = dndx(n); !  xy
      end do

      ! Compute components of the average stress and strain tensors in each pixel
      ! Initialize global accumulators
      strxx = 0.d0; stryy = 0.d0; strzz = 0.d0
      strxz = 0.d0; stryz = 0.d0; strxy = 0.d0
      sxx = 0.d0;   syy = 0.d0;   szz = 0.d0
      sxz = 0.d0;   syz = 0.d0;   sxy = 0.d0
	  
	  ! ======================================================================
      ! ---- 2. Parallel Loop over Pixels (Stress/Strain Calculation) --------
      ! Gather nodal displacements, apply periodic BC jumps, and contract with 
      ! the stiffness tensor.
      ! ======================================================================
      !$omp parallel do collapse(3) default(shared) private(m_local,mm,n8,n3,n,ibl,uu,str11_t,str22_t,str33_t,str13_t,str23_t,str12_t, s11_t,s22_t,s33_t,s13_t,s23_t,s12_t, pixm) reduction(+:strxx,stryy,strzz,strxz,stryz,strxy,sxx,syy,szz,sxz,syz,sxy) schedule(guided)
      do k_local = 1, nz
        do j_local = 1, ny
          do i_local = 1, nx
		  
            m_local = (k_local-1)*nxy + (j_local-1)*nx + i_local

            ! Gather Connectivity (Mapping 27 potential neighbors to 8 element nodes)
            ibl( 1) = e3d%ib(m_local,  1)
			ibl( 2) = e3d%ib(m_local,  2)
            ibl( 3) = e3d%ib(m_local,  3)
            ibl(17) = e3d%ib(m_local, 17)
			ibl(18) = e3d%ib(m_local, 18)
            ibl(19) = e3d%ib(m_local, 19)
			ibl(26) = e3d%ib(m_local, 26)

            ! Load local nodal displacements
            uu = 0.d0
            do mm = 1, ndof
              uu(1,mm) = e3d%u(m_local,mm)
              uu(2,mm) = e3d%u(ibl(3 ),mm)
              uu(3,mm) = e3d%u(ibl(2 ),mm)
              uu(4,mm) = e3d%u(ibl(1 ),mm)
              uu(5,mm) = e3d%u(ibl(26),mm)
              uu(6,mm) = e3d%u(ibl(19),mm)
              uu(7,mm) = e3d%u(ibl(18),mm)
              uu(8,mm) = e3d%u(ibl(17),mm)
            end do

            ! Apply Periodic Boundary Conditions
            if(i_local == nx) then
              uu(2,1) = uu(2,1) + exx * dble(nx)
              uu(3,1) = uu(3,1) + exx * dble(nx)
              uu(6,1) = uu(6,1) + exx * dble(nx)
              uu(7,1) = uu(7,1) + exx * dble(nx)
              uu(2,2) = uu(2,2) + exy * dble(nx)
              uu(3,2) = uu(3,2) + exy * dble(nx)
              uu(6,2) = uu(6,2) + exy * dble(nx)
              uu(7,2) = uu(7,2) + exy * dble(nx)
              uu(2,3) = uu(2,3) + exz * dble(nx)
              uu(3,3) = uu(3,3) + exz * dble(nx)
              uu(6,3) = uu(6,3) + exz * dble(nx)
              uu(7,3) = uu(7,3) + exz * dble(nx)
            end if
            if(j_local == ny) then
              uu(3,1) = uu(3,1) + exy * dble(ny)
              uu(4,1) = uu(4,1) + exy * dble(ny)
              uu(7,1) = uu(7,1) + exy * dble(ny)
              uu(8,1) = uu(8,1) + exy * dble(ny)
              uu(3,2) = uu(3,2) + eyy * dble(ny)
              uu(4,2) = uu(4,2) + eyy * dble(ny)
              uu(7,2) = uu(7,2) + eyy * dble(ny)
              uu(8,2) = uu(8,2) + eyy * dble(ny)
              uu(3,3) = uu(3,3) + eyz * dble(ny)
              uu(4,3) = uu(4,3) + eyz * dble(ny)
              uu(7,3) = uu(7,3) + eyz * dble(ny)
              uu(8,3) = uu(8,3) + eyz * dble(ny)
            end if
            if(k_local == nz) then
              uu(5,1) = uu(5,1) + exz * dble(nz)
              uu(6,1) = uu(6,1) + exz * dble(nz)
              uu(7,1) = uu(7,1) + exz * dble(nz)
              uu(8,1) = uu(8,1) + exz * dble(nz)
              uu(5,2) = uu(5,2) + eyz * dble(nz)
              uu(6,2) = uu(6,2) + eyz * dble(nz)
              uu(7,2) = uu(7,2) + eyz * dble(nz)
              uu(8,2) = uu(8,2) + eyz * dble(nz)
              uu(5,3) = uu(5,3) + ezz * dble(nz)
              uu(6,3) = uu(6,3) + ezz * dble(nz)
              uu(7,3) = uu(7,3) + ezz * dble(nz)
              uu(8,3) = uu(8,3) + ezz * dble(nz)
            end if

			! Reset local thread accumulators
            str11_t = 0.d0; str22_t = 0.d0; str33_t = 0.d0
            str13_t = 0.d0; str23_t = 0.d0; str12_t = 0.d0
            s11_t = 0.d0;   s22_t = 0.d0;   s33_t = 0.d0
            s13_t = 0.d0;   s23_t = 0.d0;   s12_t = 0.d0

            pixm = e3d%pix(m_local)
			
			! --- Contraction: Stress = C * B * u ---
            !$omp simd
            do n3 = 1, ndof
              do n8 = 1, nnode_fe
                s11_t = s11_t + es(1,n8,n3)*uu(n8,n3)
                s22_t = s22_t + es(2,n8,n3)*uu(n8,n3)
                s33_t = s33_t + es(3,n8,n3)*uu(n8,n3)
                s13_t = s13_t + es(4,n8,n3)*uu(n8,n3)
                s23_t = s23_t + es(5,n8,n3)*uu(n8,n3)
                s12_t = s12_t + es(6,n8,n3)*uu(n8,n3)
                do n = 1, 6
                  str11_t = str11_t + e3d%cmod(pixm,1,n)*es(n,n8,n3)*uu(n8,n3)
                  str22_t = str22_t + e3d%cmod(pixm,2,n)*es(n,n8,n3)*uu(n8,n3)
                  str33_t = str33_t + e3d%cmod(pixm,3,n)*es(n,n8,n3)*uu(n8,n3)
                  str13_t = str13_t + e3d%cmod(pixm,4,n)*es(n,n8,n3)*uu(n8,n3)
                  str23_t = str23_t + e3d%cmod(pixm,5,n)*es(n,n8,n3)*uu(n8,n3)
                  str12_t = str12_t + e3d%cmod(pixm,6,n)*es(n,n8,n3)*uu(n8,n3)
                end do
              end do
            end do

            ! Accumulate thread results to global sums (Reduction)
            strxx = strxx + str11_t
            stryy = stryy + str22_t
            strzz = strzz + str33_t
            strxz = strxz + str13_t
            stryz = stryz + str23_t
            strxy = strxy + str12_t
            sxx   = sxx   + s11_t
            syy   = syy   + s22_t
            szz   = szz   + s33_t
            sxz   = sxz   + s13_t
            syz   = syz   + s23_t
            sxy   = sxy   + s12_t

          end do
        end do
      end do
      !$omp end parallel do

      ! ======================================================================
      ! ---- 3. Volume Averaging and Persistence -----------------------------
      ! Normalize sums by the total number of samples (voxels) and store back
      ! in the global data structure.
      ! ======================================================================
	  ! Volume average
      strxx = strxx / dble(ns)
      stryy = stryy / dble(ns)
      strzz = strzz / dble(ns)
      strxz = strxz / dble(ns)
      stryz = stryz / dble(ns)
      strxy = strxy / dble(ns)
	  
      sxx = sxx / dble(ns)
      syy = syy / dble(ns)
      szz = szz / dble(ns)
      sxz = sxz / dble(ns)
      syz = syz / dble(ns)
      sxy = sxy / dble(ns)

      ! Assign back to struct
      e3d%strxx = strxx
      e3d%stryy = stryy
      e3d%strzz = strzz
      e3d%strxz = strxz
      e3d%stryz = stryz
      e3d%strxy = strxy
	  
      e3d%sxx   = sxx
      e3d%syy   = syy
      e3d%szz   = szz
      e3d%sxz   = sxz
      e3d%syz   = syz
      e3d%sxy   = sxy
      
      end subroutine stress_OpenMP
!-------------------------------------------------------------------------     
      subroutine assig(e3d, prob)
	  ! ==============================================================================
      ! SUBROUTINE: assig
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Calculates the volume fraction (probability) of each material phase present 
      !   in the microstructure. It scans the entire voxel grid and normalizes 
      !   phase counts against the total number of samples.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(in)): Derived type containing the voxel phase map (pix).
      !   - prob (intent(out)): An array of size 'nphase' populated with the 
      !     calculated volume fractions [0.0 to 1.0].
      !
      ! 'e3d' COMPONENTS MODIFIED/ACCESSED:
      !   - e3d%pix (intent(in)): The phase identifier for every voxel in the grid.
      !
      ! METHODOLOGY:
      !   - Discrete Counting: Iterates through the global sample set (ns) and 
      !     increments a counter corresponding to the phase ID at each voxel.
      !   - Validation: Includes a bounds check to ensure voxel IDs fall within 
      !     the defined phase range [1, nphase].
      !   - Normalization: Divides the final counts by the total number of voxels 
      !     to transform absolute counts into relative volume fractions.
      ! ==============================================================================	  
      use elas3d_mod
      implicit none
      
      ! Arguments
      type(elas3d_data_type), intent(in) :: e3d
      real(dp), intent(out) :: prob(nphase)
      
      ! Locals
      integer :: m

      ! ======================================================================
      ! ---- 1. Initialization -----------------------------------------------
      ! Ensure the output array is zeroed before starting the summation.
      ! ======================================================================
      prob(:) = 0.d0

      ! ======================================================================
      ! ---- 2. Phase Counting -----------------------------------------------
      ! Traverse the 1D phase map array and accumulate occurrences.
      ! ======================================================================
      do m = 1, ns
        if(e3d%pix(m) >= 1 .and. e3d%pix(m) <= nphase) then
          prob(e3d%pix(m)) = prob(e3d%pix(m)) + 1.d0
        end if
      end do

      ! ======================================================================
      ! ---- 3. Fractional Normalization -------------------------------------
      ! Scale the counts by the total volume (ns) to get volume fractions.
      ! ======================================================================
      prob(:) = prob(:) / dble(ns)


      end subroutine assig	  
!-------------------------------------------------------------------------  
      subroutine rotatestiffnessby(e3d)
	  ! ==============================================================================
      ! SUBROUTINE: rotatestiffnessby
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Computes the rotated 6x6 stiffness tensors for each polycrystalline phase.
      !   It transforms the intrinsic crystal-frame stiffness into the global 
      !   laboratory frame using Rodrigues orientation vectors and Bond 
      !   transformation matrices.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(inout)): Derived type containing orientation data and 
      !     the destination array for rotated tensors.
      !
      ! 'e3d' COMPONENTS MODIFIED/ACCESSED:
      !   - e3d%orientation (intent(in)): Rodrigues vectors for each phase.
      !   - e3d%rotatedstiffness (intent(out)): Storage for the 6x6 rotated tensors.
      !
      ! METHODOLOGY:
      !   - Rotation Logic: Uses the Rodrigues formula to build a 3x3 rotation 
      !     matrix (Q) from the orientation vector (r).
      !   - Bond Transformation: Constructs a 6x6 "newr" matrix (T) to handle the 
      !     transformation of 4th-order tensors in Voigt notation.
      !   - Math: Performs the triple matrix product: $C_{global} = T^T C_{crystal} T$.
      !   - Note: Includes proper engineering shear strain factors (2.0) in the 
      !     lower-left quadrant of the Bond matrix to ensure energetic consistency.
      ! ==============================================================================
      use elas3d_mod
      implicit none

      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d
      
      ! Locals
      integer :: i, j, k, l, g
      real(dp) :: term1, term2, r(3), rdotr
	  real(dp) :: StiffnessInCrystalFrame(6,6)
	  real(dp) :: OrientationMatrix(3,3)
      real(dp) :: newr(6,6)
      real(dp) :: rotated(6,6)
	  real(dp) :: Q(3,3)

	  
	  ! ======================================================================
      ! ---- 1. Initialize Crystal Frame Stiffness ---------------------------
      ! Define the base orthotropic/cubic stiffness for the material.
      ! ======================================================================
      StiffnessInCrystalFrame(:,:) = 0.d0
      StiffnessInCrystalFrame(1,1) = C11_local
      StiffnessInCrystalFrame(2,2) = C11_local
      StiffnessInCrystalFrame(3,3) = C11_local
      StiffnessInCrystalFrame(1,2) = C12_local
      StiffnessInCrystalFrame(1,3) = C12_local
      StiffnessInCrystalFrame(2,3) = C12_local
      StiffnessInCrystalFrame(2,1) = C12_local
      StiffnessInCrystalFrame(3,1) = C12_local
      StiffnessInCrystalFrame(3,2) = C12_local
      StiffnessInCrystalFrame(4,4) = C44_local
      StiffnessInCrystalFrame(5,5) = C44_local
      StiffnessInCrystalFrame(6,6) = C44_local

      ! ======================================================================
      ! ---- 2. Loop Over Phases (Orientation Transformation) ----------------
      ! Iterate through each phase to apply its specific orientation vector.
      ! ======================================================================
      do g = 1, nphase-1
	  
	    ! Load Orientation
        r(1:3) = e3d%orientation((g-1)*3+1:(g-1)*3+3)    
		              
        ! Build 3x3 Orientation Matrix via Rodrigues Formula
        rdotr = sum(r(:)**2)
        term1 = 1.d0 - rdotr
        term2 = 1.d0 + rdotr
		
        OrientationMatrix(:,:) = 0.d0
        OrientationMatrix(1,1) = term1
        OrientationMatrix(2,2) = term1
        OrientationMatrix(3,3) = term1
        
		! Rodriques Formula
        OrientationMatrix = OrientationMatrix + 2.d0 * spread(r, 2, 3) * spread(r, 1, 3)
       
	    ! Off-diagonal terms (incorporating cross product components)
        OrientationMatrix(1,2) = OrientationMatrix(1,2) - 2.d0*r(3)
        OrientationMatrix(1,3) = OrientationMatrix(1,3) + 2.d0*r(2)
        OrientationMatrix(2,3) = OrientationMatrix(2,3) - 2.d0*r(1)
        OrientationMatrix(2,1) = OrientationMatrix(2,1) + 2.d0*r(3)
        OrientationMatrix(3,1) = OrientationMatrix(3,1) - 2.d0*r(2)
        OrientationMatrix(3,2) = OrientationMatrix(3,2) + 2.d0*r(1)
        
        ! Normalize Rotation Matrix (Handle Identity for zero rotation)
        if (term2 > 1.0d-12) then
           term2 = 1.0d0 / term2
           OrientationMatrix(:,:) = OrientationMatrix(:,:) * term2
        else
           OrientationMatrix = 0.d0
           OrientationMatrix(1,1) = 1.0d0
           OrientationMatrix(2,2) = 1.0d0
           OrientationMatrix(3,3) = 1.0d0
        end if
		
        Q = OrientationMatrix
        ! ======================================================================
        ! ---- 3. Build 6x6 Bond Transformation Matrix (Strain Basis) ----------
        ! Constructs 'newr' to transform Voigt-notated tensors.
        ! ======================================================================
        ! Row 1		
        newr(1,1) = Q(1,1)**2
        newr(1,2) = Q(1,2)**2
        newr(1,3) = Q(1,3)**2
        newr(1,4) = Q(1,2)*Q(1,3)
        newr(1,5) = Q(1,3)*Q(1,1)
        newr(1,6) = Q(1,1)*Q(1,2)
 	    ! Row 2   
        newr(2,1) = Q(2,1)**2
        newr(2,2) = Q(2,2)**2
        newr(2,3) = Q(2,3)**2
        newr(2,4) = Q(2,2)*Q(2,3)
        newr(2,5) = Q(2,3)*Q(2,1)
        newr(2,6) = Q(2,1)*Q(2,2)
 	    ! Row 3   
        newr(3,1) = Q(3,1)**2
        newr(3,2) = Q(3,2)**2
        newr(3,3) = Q(3,3)**2
        newr(3,4) = Q(3,2)*Q(3,3)
        newr(3,5) = Q(3,3)*Q(3,1)
        newr(3,6) = Q(3,1)*Q(3,2)
    
        ! Rows 4 (With Engineering Strain Multipliers)	
        newr(4,1) = 2.d0 * Q(3,1)*Q(1,1)
        newr(4,2) = 2.d0 * Q(3,2)*Q(1,2)
        newr(4,3) = 2.d0 * Q(3,3)*Q(1,3)
        newr(4,4) = Q(1,2)*Q(3,3) + Q(1,3)*Q(3,2)
        newr(4,5) = Q(1,3)*Q(3,1) + Q(1,1)*Q(3,3)
        newr(4,6) = Q(1,1)*Q(3,2) + Q(1,2)*Q(3,1)
    
	    ! Rows 5
        newr(5,1) = 2.d0 * Q(2,1)*Q(3,1)
        newr(5,2) = 2.d0 * Q(2,2)*Q(3,2)
        newr(5,3) = 2.d0 * Q(2,3)*Q(3,3)
        newr(5,4) = Q(2,2)*Q(3,3) + Q(2,3)*Q(3,2)
        newr(5,5) = Q(2,1)*Q(3,3) + Q(2,3)*Q(3,1)
        newr(5,6) = Q(2,2)*Q(3,1) + Q(2,1)*Q(3,2)	
    
	    ! Rows 6
        newr(6,1) = 2.d0 * Q(1,1)*Q(2,1)
        newr(6,2) = 2.d0 * Q(1,2)*Q(2,2)
        newr(6,3) = 2.d0 * Q(1,3)*Q(2,3)
        newr(6,4) = Q(1,2)*Q(2,3) + Q(1,3)*Q(2,2)
        newr(6,5) = Q(1,3)*Q(2,1) + Q(1,1)*Q(2,3)
        newr(6,6) = Q(1,1)*Q(2,2) + Q(1,2)*Q(2,1)      

        ! ======================================================================
        ! ---- 4. Execute Rotation and Global Storage --------------------------
        ! $C_{global} = T^T \cdot C_{crystal} \cdot T$
        ! ======================================================================
        rotated = matmul(matmul(transpose(newr), StiffnessInCrystalFrame), newr)	
		
		do i = 1, 6
		  do l = 1, 6
            e3d%rotatedstiffness((g-1)*36 + (i-1)*6 + l) = rotated(i, l)
          end do
        end do
		
      end do
      
      end subroutine rotatestiffnessby	  
!-------------------------------------------------------------------------
      subroutine stress_fullfield_OpenMP(e3d) 
	  ! ==============================================================================
      ! SUBROUTINE: stress_fullfield_OpenMP
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Computes the full-field Cauchy stress tensor and von Mises equivalent 
      !   stress for every voxel in the microstructure. Results are exported to 
      !   an HDF5 file for advanced post-processing and visualization.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(inout)): Derived type containing displacements, phase maps,
      !     and the storage arrays for the full-field results.
      !
      ! 'e3d' COMPONENTS MODIFIED/ACCESSED:
      !   - e3d%u (intent(in)): Nodal displacement field.
      !   - e3d%stress_field (intent(out)): Local stress tensor (ns x 6).
      !   - e3d%vm (intent(out)): von Mises equivalent stress field.
      !
      ! METHODOLOGY:
      !   - Voxel-wise Calculation: Similar to 'stress_OpenMP', but stores results
      !     locally for every element instead of just volume averaging.
      !   - von Mises: Computes the second invariant of the deviatoric stress tensor.
      !   - I/O: Utilizes the HDF5 library for high-performance binary data export,
      !     storing phase maps, von Mises stress, and the full stress tensor.
      ! ==============================================================================
      use elas3d_mod
      use hdf5
      implicit none
  
      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d

      ! Loop variables
      integer :: i, j, k, m, nxy, mm, n8, n3, n
	  integer :: m_l, k_l, j_l, i_l 
	  
      ! Element Strain-Displacement Matrix (6 components, 8 nodes, 3 DOFs)	  
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
      real(dp) :: es(6,nnode_fe,ndof)         

      ! Thread-private arrays
      real(dp) :: uu(nnode_fe,ndof)
      real(dp) :: eps(6)	
      integer :: ibl(27)	  

      ! Macroscopic strains 
      real(dp) :: exx, eyy, ezz, exz, eyz, exy  

      ! HDF5 variables
      integer(HID_T)    :: file_id, dset_id, space_id
      integer(HSIZE_T)  :: dims(2)
      integer           :: hdferr   

      ! ======================================================================
      ! ---- 1. Initialization and B-Matrix Setup ----------------------------
      ! ======================================================================
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy  
  
      ! Physical sizes
      nxy = nx * ny  

      ! Trilinear shape function derivatives at voxel center
      dndx = [-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0,-0.25d0]
      dndy = [-0.25d0,-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0]
      dndz = [-0.25d0,-0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0]
      
      ! Build Strain-Displacement Matrix (6 components x 8 nodes x 3 DOFs)
      ! Mapping: 1=xx, 2=yy, 3=zz, 4=xz, 5=yz, 6=xy
      es = 0.d0
      do n = 1,8
        es(1,n,1) = dndx(n)
        es(2,n,2) = dndy(n)
        es(3,n,3) = dndz(n)
        es(4,n,1) = dndz(n);    es(4,n,3) = dndx(n)
        es(5,n,2) = dndz(n);    es(5,n,3) = dndy(n)
        es(6,n,1) = dndy(n);    es(6,n,2) = dndx(n)
      end do


      ! ======================================================================
      ! ---- 2. Parallel Full-Field Stress Calculation -----------------------
      ! ======================================================================
      ! Main per-voxel calculation 
      !$omp parallel do collapse(3) schedule(guided) default(shared) private(k_l,j_l,i_l,m_l,mm,n8,n3,ibl,uu,eps) 
      do k_l = 1, nz
	    do j_l = 1, ny
		  do i_l = 1, nx
        
		    m_l = (k_l-1)*nxy + (j_l-1)*nx + i_l


            ! Gather Connectivity (Mapping 27 potential neighbors to 8 element nodes)
            ibl( 1) = e3d%ib(m_l,  1); 
		    ibl( 2) = e3d%ib(m_l,  2); 
		    ibl( 3) = e3d%ib(m_l,  3); 
            ibl(17) = e3d%ib(m_l, 17); 
		    ibl(18) = e3d%ib(m_l, 18); 
		    ibl(19) = e3d%ib(m_l, 19); 
		    ibl(26) = e3d%ib(m_l, 26);
        
            ! Load local nodal displacements
            uu = 0.d0
            do mm = 1, ndof
              uu(1,mm) = e3d%u(m_l,    mm)
              uu(2,mm) = e3d%u(ibl(3 ),mm)
              uu(3,mm) = e3d%u(ibl(2 ),mm)
              uu(4,mm) = e3d%u(ibl(1 ),mm)
              uu(5,mm) = e3d%u(ibl(26),mm)
              uu(6,mm) = e3d%u(ibl(19),mm)
              uu(7,mm) = e3d%u(ibl(18),mm)
              uu(8,mm) = e3d%u(ibl(17),mm)
            end do
    
	
            ! Apply Periodic Boundary Conditions
            if (i_l == nx) then
              uu(2,1) = uu(2,1) + exx * dble(nx)
              uu(3,1) = uu(3,1) + exx * dble(nx)
              uu(6,1) = uu(6,1) + exx * dble(nx)
              uu(7,1) = uu(7,1) + exx * dble(nx)
              uu(2,2) = uu(2,2) + exy * dble(nx)
              uu(3,2) = uu(3,2) + exy * dble(nx)
              uu(6,2) = uu(6,2) + exy * dble(nx)
              uu(7,2) = uu(7,2) + exy * dble(nx)
              uu(2,3) = uu(2,3) + exz * dble(nx)
              uu(3,3) = uu(3,3) + exz * dble(nx)
              uu(6,3) = uu(6,3) + exz * dble(nx)
              uu(7,3) = uu(7,3) + exz * dble(nx)
            end if
            if (j_l == ny) then
              uu(3,1) = uu(3,1) + exy * dble(ny)
              uu(4,1) = uu(4,1) + exy * dble(ny)
              uu(7,1) = uu(7,1) + exy * dble(ny)
              uu(8,1) = uu(8,1) + exy * dble(ny)
              uu(3,2) = uu(3,2) + eyy * dble(ny)
              uu(4,2) = uu(4,2) + eyy * dble(ny)
              uu(7,2) = uu(7,2) + eyy * dble(ny)
              uu(8,2) = uu(8,2) + eyy * dble(ny)
              uu(3,3) = uu(3,3) + eyz * dble(ny)
              uu(4,3) = uu(4,3) + eyz * dble(ny)
              uu(7,3) = uu(7,3) + eyz * dble(ny)
              uu(8,3) = uu(8,3) + eyz * dble(ny)
            end if
            if (k_l == nz) then
              uu(5,1) = uu(5,1) + exz * dble(nz)
              uu(6,1) = uu(6,1) + exz * dble(nz)
              uu(7,1) = uu(7,1) + exz * dble(nz)
              uu(8,1) = uu(8,1) + exz * dble(nz)
              uu(5,2) = uu(5,2) + eyz * dble(nz)
              uu(6,2) = uu(6,2) + eyz * dble(nz)
              uu(7,2) = uu(7,2) + eyz * dble(nz)
              uu(8,2) = uu(8,2) + eyz * dble(nz)
              uu(5,3) = uu(5,3) + ezz * dble(nz)
              uu(6,3) = uu(6,3) + ezz * dble(nz)
              uu(7,3) = uu(7,3) + ezz * dble(nz)
              uu(8,3) = uu(8,3) + ezz * dble(nz)
            end if         


            ! Compute strain vector (eps) and local stress via Hooke's Law
            eps = 0.d0
            do n3 = 1, ndof
              do n8 = 1, nnode_fe
                eps(1) = eps(1) + es(1,n8,n3)*uu(n8,n3)
                eps(2) = eps(2) + es(2,n8,n3)*uu(n8,n3)
                eps(3) = eps(3) + es(3,n8,n3)*uu(n8,n3)
                eps(4) = eps(4) + es(4,n8,n3)*uu(n8,n3)
                eps(5) = eps(5) + es(5,n8,n3)*uu(n8,n3)
                eps(6) = eps(6) + es(6,n8,n3)*uu(n8,n3)
              end do
            end do

            ! Store the Cauchy stress tensor for this voxel
            e3d%stress_field(m_l,:) = matmul(e3d%cmod(e3d%pix(m_l),1:6,1:6), eps)
      
          end do
		end do
      end do
      !$omp end parallel do

      
      ! ======================================================================
      ! ---- 3. von Mises Equivalent Stress Calculation ----------------------
      ! ======================================================================  
      !$omp parallel do schedule(guided)
      do m = 1, ns
        e3d%vm(m) = sqrt(0.5d0*((e3d%stress_field(m,1)-e3d%stress_field(m,2))**2 + &
                                (e3d%stress_field(m,2)-e3d%stress_field(m,3))**2 + &
                                (e3d%stress_field(m,3)-e3d%stress_field(m,1))**2 + &
                         6.0d0*(e3d%stress_field(m,4)**2 + &
                                e3d%stress_field(m,5)**2 + &
						        e3d%stress_field(m,6)**2)))
      end do
      !$omp end parallel do


      ! ======================================================================
      ! ---- 4. HDF5 Data Export ---------------------------------------------
      ! Writes voxel data, von Mises stress, and full tensors to a binary file.
      ! ======================================================================  
      call h5open_f(hdferr)
      ! Create file
      call h5fcreate_f('fullfield_poly.h5', H5F_ACC_TRUNC_F, file_id, hdferr)
	  
      ! Create dataspace
      dims(1) = ns	  
      call h5screate_simple_f(1, dims, space_id, hdferr)
	  
	  ! A. Create dataset (integers), pix(m)
      call h5dcreate_f(file_id, "pix", H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
	  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, e3d%pix, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
	  
	  ! B. Create dataset (doubles), vm(m)
      call h5dcreate_f(file_id, "vm", H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, e3d%vm, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
	  ! Cleanup
      call h5sclose_f(space_id, hdferr)
	  
      ! Write Full Stress Tensor (ns x 6)
      ! C. Create dataset for stress tensor (ns x 6)
      dims(1) = ns
      dims(2) = 6
      call h5screate_simple_f(2, dims, space_id, hdferr)

      call h5dcreate_f(file_id, "stress_tensor", H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, e3d%stress_field, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
	  
	  ! Cleanup
	  call h5sclose_f(space_id, hdferr)	  	  
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

      
      end subroutine stress_fullfield_OpenMP    
!-------------------------------------------------------------------------    
      subroutine ppixel_hdf5(e3d)
	  ! ==============================================================================
      ! SUBROUTINE: ppixel_hdf5
      ! ------------------------------------------------------------------------------
      ! PURPOSE:
      !   Reads the input microstructure (phase map) and crystallographic 
      !   orientation vectors from an external HDF5 file ('input_structure_poly.h5'). 
      !   This serves as the initialization step for the digital RVE.
      !
      ! DUMMY ARGUMENTS:
      !   - e3d (intent(inout)): Derived type where the 'pix' (phase IDs) and 
      !     'orientation' (Rodrigues vectors) arrays will be populated.
      !
      ! 'e3d' COMPONENTS MODIFIED/ACCESSED:
      !   - e3d%pix (intent(out)): Populated with phase integers for each voxel.
      !   - e3d%orientation (intent(out)): Populated with Rodrigues-Frank vectors.
      !
      ! METHODOLOGY:
      !   - File I/O: Uses the HDF5 Fortran interface to open and read datasets.
      !   - Verification: Performs strict rank and dimension checking to ensure 
      !     the input file matches the current simulation's grid size (ns) and 
      !     phase count (nphase).
      !   - Boundary Check: Validates that all imported phase labels fall within 
      !     the physical range [1, nphase].
      ! ==============================================================================
      use elas3d_mod
      use hdf5
      implicit none

      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d

      ! Locals
      integer(HID_T)   :: file_id, dset_id, space_id
      integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
      integer           :: rank, hdferr
	  integer(HSIZE_T)  :: total_elems
      integer(HSIZE_T)  :: odims(1)

	  
	  ! ======================================================================
      ! ---- 1. Open HDF5 Library and File -----------------------------------
      ! ======================================================================
      call h5open_f(hdferr)
      call h5fopen_f('input_structure_poly.h5', H5F_ACC_RDONLY_F, file_id, hdferr)
	  if (hdferr /= 0) stop 'Error opening HDF5 file'
	  
	  ! ======================================================================
      ! ---- 2. Read 'pix' Dataset (Microstructure Phase Map) ----------------
      ! ======================================================================
      call h5dopen_f(file_id, 'pix', dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)

      
      ! Query and verify dataset dimensions
      call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
      allocate(dims(rank), maxdims(rank))
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)

	  total_elems = product(dims) 
      if (total_elems /= ns) then
        print *, 'ERROR: HDF5 pix size (', total_elems, ') does not match simulation ns (', ns, ')'
        stop
      end if

      ! Transfer data to global structure
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, e3d%pix, dims, hdferr)
	  
	  ! Cleanup pix-related handles
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      deallocate(dims, maxdims)
	  
	  ! ======================================================================
      ! ---- 3. Read 'orientation' Dataset (Rodrigues Vectors) ---------------
      ! ======================================================================
      call h5dopen_f(file_id, 'orientation', dset_id, hdferr)
      ! Get rank and dims for orientation
      call h5dget_space_f(dset_id, space_id, hdferr)
	  
	  ! Get dimensions to verify size (rank)
      call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
      allocate(dims(rank), maxdims(rank))
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
	  	  
      ! Verify size matches allocated orientation array
	  total_elems = product(dims)
      if (total_elems /= size(e3d%orientation)) then
          print *, 'ERROR: HDF5 orientation size (', total_elems, ') mismatch with allocated (', size(e3d%orientation), ')'
          stop
      end if
	  
	  ! Transfer data using 1D dimension specifier
      odims(1) = size(e3d%orientation)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, e3d%orientation, odims, hdferr)
	  
	  ! Cleanup 'orientation' handles
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      deallocate(dims, maxdims)


      ! Finished
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

      ! Ensure phase labels are consistent with the simulation parameters.
      if (minval(e3d%pix) < 1 .or. maxval(e3d%pix) > nphase) then
        print *, "ERROR: Phase label in pix out of bounds"
        stop
      end if

      end subroutine ppixel_hdf5      