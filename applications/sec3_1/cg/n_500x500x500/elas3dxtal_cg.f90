!------------------------------------------------------------------------------
! ELAS3D-Xtal: Conjugate Gradient Solver for 3D Crystal Elasticity
! ------------------------------------------------------------------------------
! This program solves the 3D linear elastic response of a voxelized multi-phase
! microstructure (e.g., XCT-derived defect-containing metal) using a regular 
! grid finite element method. The code supports both isotropic and anisotropic 
! phases and applies periodic boundary conditions. The stress and strain
! response to user-prescribed average strains is computed via a conjugate
! gradient energy minimization algorithm, heavily accelerated using OpenMP.
!    
! This program is based on the original work:
!   Garboczi, E. J. (1998). Finite element and finite difference programs for 
!   computing the linear electric and elastic properties of digital images of 
!   random materials. NISTIR 6269. National Institute of Standards and Technology.    
!
! Major capabilities:
!   - Reads microstructure (phase) and orientation data via HDF5.
!   - Assembles periodic 3D finite element stiffness matrices.
!   - Supports isotropic matrix/inclusion or polycrystalline matrix/inclusion.
!   - Applies macroscopic loading (strain-controlled).
!   - Solves for equilibrium using OpenMP-parallelized conjugate gradient minimization.
!   - Extracts and writes full-field stress and von Mises fields for all voxels.
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
! Dates & Revisions:
!   July 28, 2025  - Initial parallel formulation
!   Aug 22, 2025   - Memory Optimizations, HDF5 Integration
!   Aug 28, 2025   - OpenMP chunking iteration implemented (dembx_OpenMP)
! ------------------------------------------------------------------------------
      module elas3d_mod
       use, intrinsic :: iso_fortran_env, dp=>real64 
       implicit none
       
       ! ---- Domain/Grid/Material/Parameters as parameters ----
       integer, parameter :: md = 500   
       integer, parameter :: nx = md     
       integer, parameter :: ny = md     
       integer, parameter :: nz = md     
       integer, parameter :: ns = nx * ny * nz   

       ! ---- Microstructure ----
       integer, parameter :: n_grains = 1        
       integer, parameter :: nphase = n_grains+1 
       integer, parameter :: nphmax = nphase     
       integer, parameter :: nfaces = 27         
       integer, parameter :: ndof = 3            
	   integer, parameter :: nnode_fe = 8        
       integer, parameter :: ngauss   = 3       

       ! ---- Components of Elastic stiffness tensor (FCC) ----
       real(dp), parameter :: C11_local = 256.75d9   ! (Pa)
       real(dp), parameter :: C12_local = 106.59d9   ! (Pa)
       real(dp), parameter :: C44_local = 75.08d9    ! (Pa)
             
       ! ---- Isotropic Material Constnats ----
	   real(dp), parameter :: K0 = (C11_local + 2.0d0*C12_local)/3.0d0
       real(dp), parameter :: G0 = (C11_local - C12_local + 3.0d0*C44_local)/5.0d0	   
       real(dp), parameter :: E0 = 9.0d0*K0*G0/(3.0d0*K0 + G0)
       real(dp), parameter :: nu0 = (3.0d0*K0 - 2.0d0*G0) / (2.0d0*(3.0d0*K0+G0))  ! Matrix
       real(dp), parameter :: E1 = 0.0d0        ! Inclusion
       real(dp), parameter :: nu1 = 0.0d0       ! Inclusion        
       integer, parameter :: flag_m = 0         ! 0: isotropic matrix, 1: polycrystal	   
              
       ! ---- Applied macroscopic loading ----
       real(dp), parameter :: aml     =  1.0e-3
       real(dp), parameter :: aml_exx =  0.00d0
       real(dp), parameter :: aml_eyy =  0.00d0       
       real(dp), parameter :: aml_ezz =  aml
       real(dp), parameter :: aml_exz =  0.00d0
       real(dp), parameter :: aml_eyz =  0.00d0
       real(dp), parameter :: aml_exy =  0.00d0       
       
       ! ---- Iterations ----
       integer, parameter :: kmax = 1              
       integer, parameter :: ldemb = 10000         
       integer, parameter :: n_iter = kmax * ldemb 
       integer, parameter :: block_size = 8192        
       
       ! ---- Convergence ----
       real, parameter :: tol = 1.0d-8

       ! ---- Derived Type: the structure ----
       type :: elas3d_data_type
         ! FE/solver variables
         real(dp), allocatable  :: u(:,:), gb(:,:), b(:,:)
		 real(dp), allocatable  :: h(:,:)
         real(dp), allocatable  :: cmod(:,:,:), dk(:,:,:,:,:)
         integer, allocatable   :: ib(:,:), pix(:)
         real(dp) :: exx, eyy, ezz, exz, eyz, exy
         real(dp) :: C
         real(dp) :: strxx, stryy, strzz, strxz, stryz, strxy
         real(dp) :: sxx, syy, szz, sxz, syz, sxy
       
         ! Per-voxel full field output
         real(dp), allocatable  :: vm(:), stress_field(:,:)
       
         ! Microstructural/iteration variables
         real(dp), allocatable  :: orientation(:)
         real(dp), allocatable  :: rotatedstiffness(:)
         real(dp), allocatable  :: energies(:), ggs(:)
         integer, allocatable  :: cgiter(:)

       end type elas3d_data_type
       
      end module elas3d_mod
    
!=========================================================================    
      program elas3d
      use elas3d_mod
      use omp_lib   
      implicit none

      ! Variable declarations   
      integer :: i1, j1, k1, m1
      integer :: i, j, k, m, n
      integer :: micro, npoints, ltot, kkk, Lstep
      integer :: ios, nxy
      integer :: im(nfaces), jm(nfaces), km(nfaces)
      real(dp) :: gg, utot, youngs
      real(dp) :: x, y, z, elapsed_time, gg0, phmod(2,2)
      real(dp), allocatable :: phasemod(:,:), prob(:)
      integer :: t1, t2, tc 
      integer :: lhist
      
      type(elas3d_data_type) :: e3d
      
  
      ! --- ALLOCATE FIELDS INSIDE STRUCT ---    
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

      
      ! INITIALIZE FIELDS
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
      
      ! FILE OPERATIONS 
      open(unit=7, file='output_iso.txt', status='replace', iostat=ios)
      
      ! ---- Timing Start (Fortran: use `system_clock`) ----
      call system_clock(t1, tc)  
      
      ! ---- Output grid info ----
      write(7, '(a,i4,a,i4,a,i4,a,i12)') 'nx = ', nx, ', ny =', ny, ', nz =', nz, ', ns =', ns
      print '(a,i4,a,i4,a,i4,a,i12)', 'nx = ', nx, ', ny =', ny, ', nz =', nz, ', ns =', ns


      ! The parameter phasemod(i,j) is the bulk(i,1) and shear(i,2) moduli of the i'th phase.
      if (flag_m == 0) then
        do i = 1, nphase-1
           phasemod(i,1) = E0   ! matrix
           phasemod(i,2) = nu0
        end do
      end if
      phasemod(nphase,1) = E1  ! Inclusions
      phasemod(nphase,2) = nu1    
      


      ! bulk modulus (1) and shear modulus (2)
      do i = 1,nphase        
        youngs = phasemod(i,1)
        phasemod(i,1) = phasemod(i,1) / 3.d0 / (1.d0 - 2.d0*phasemod(i,2))
        phasemod(i,2) = youngs / 2.d0 / (1.d0 + phasemod(i,2))
      end do
      

      ! neighbor table, ib(m,n)
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

      ! neighbor table according to 1-d labels
      nxy = nx * ny
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            m = nxy * (k-1) + nx * (j-1) + i
            do n = 1, 27
              i1 = i + im(n)
              j1 = j + jm(n)
              k1 = k + km(n)
              if (i1 < 1)  i1 = i1 + nx
              if (i1 > nx) i1 = i1 - nx
              if (j1 < 1)  j1 = j1 + ny
              if (j1 > ny) j1 = j1 - ny
              if (k1 < 1)  k1 = k1 + nz
              if (k1 > nz) k1 = k1 - nz
              m1 = nxy * (k1-1) + nx*(j1-1) + i1
              e3d%ib(m,n) = m1
            end do
          end do
        end do
      end do
      
      
      ! Compute the average stress and strain in each microstructure.
      npoints = 1  ! Number of microstructures
      do micro = 1, npoints
        ! Read in a microstructre
        call ppixel_hdf5(e3d)
      
        print *, 'Unique phases in pix:', minval(e3d%pix), maxval(e3d%pix)
        do i = 1, nphase
          write(7, '(A,I7,A,F20.6,A,F20.6)') 'Phase ',i,' bulk=',phasemod(i,1),' shear=',phasemod(i,2)
        end do        
        
        ! Count and output the volume fractions of the different phases
        call assig(e3d, prob)
        do i = 1,nphase
          print *, 'Phase', i, 'volfrac:', prob(i)
        end do		
        do i = 1, nphase
          write(7, '(A,I7,A,F20.6)') 'Volume fraction of phase ',i,' is ',prob(i)
        end do
     

        ! Set applied strains
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
        
        ! Set up the elastic modulus variables, finite element stiffness matrices,
        ! the constant, C, and vector, b, required for computing the energy.
		if (flag_m == 0) then
		  phmod(1:2,1:2) = phasemod(1:nphase,1:2)
		else
          phmod(1:1,1:2) = phasemod(nphase:nphase,1:2)    
        end if
        call femat_OpenMP(e3d, phmod)
        deallocate(e3d%orientation, e3d%rotatedstiffness)
      
        ! Apply chosen strains as a homogeneous macroscopic strain 
        ! as the initial condition.       
        do k=1,nz; do j=1,ny; do i=1,nx
           m = nxy*(k-1)+nx*(j-1)+i
           x = dble(i-1)
           y = dble(j-1)
           z = dble(k-1)
           e3d%u(m,1) = x * e3d%exx + y * e3d%exy + z * e3d%exz
           e3d%u(m,2) = x * e3d%exy + y * e3d%eyy + z * e3d%eyz
           e3d%u(m,3) = x * e3d%exz + y * e3d%eyz + z * e3d%ezz
        end do; end do; end do
     
        
        !  RELAXATION LOOP
        ltot = 0
        lhist = 1
        
        !  Call energy to get initial energy and initial gradient
        call energy_OpenMP(e3d, utot)
        
        
        ! Initial gradient norm
        gg = sum(e3d%gb(:,:)**2)
        gg0 = gg    ! initial squared norm
        write(7,*) 'Initial energy=', utot, 'gg=', gg
        write(7,*) '   '
        call flush(7)
        print '(A,E12.5,A,E12.5)', 'Initial energy=', utot, ', Initial gg=', gg
        print *, '   '
        
        e3d%energies(lhist) = utot
        e3d%ggs(lhist)      = sqrt(gg/gg0)
        e3d%cgiter(lhist)   = 0
        
        ! Open file for writing CG history
        open(unit=99, file='cg_history_iso.txt', status='replace')
        write(99, '(I8,1X,1PE18.10,1X,1PE18.10)') e3d%cgiter(lhist), e3d%energies(lhist), e3d%ggs(lhist)
        
        do kkk=1, kmax
           ! conjugate gradient solver
           call dembx_OpenMP(e3d, Lstep, gg, gg0, kkk) 
                       

           ltot = ltot + Lstep
           lhist = lhist + 1
           
           !  Compute energy after dembx call 
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

                      
           !  If relaxation process is finished, jump out of loop
           if (sqrt(gg/gg0) < tol) exit
           
           
           !  Compute and output stresses
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

        call stress_OpenMP(e3d)


        write(7,*) 'stresses: xx, yy, zz, xz, yz, xy'
        write(7,'(6E15.6)') e3d%strxx, e3d%stryy, e3d%strzz, e3d%strxz, e3d%stryz, e3d%strxy
        write(7,*) 'strains: xx, yy, zz, xz, yz, xy'
        write(7,'(6E15.6)') e3d%sxx, e3d%syy, e3d%szz, e3d%sxz, e3d%syz, e3d%sxy

        
      end do
      
      print *, 'End of CG'
      
      !------------------------------------------
      ! TIMING STOP
      !------------------------------------------
      call system_clock(t2)
      elapsed_time = dble(t2-t1)/dble(tc)
      print *, ''
      print *, 'Calculation complete!'
      print *, 'Total calculation time: ', elapsed_time, ' seconds'        
	  write(7,*) ''
	  write(7,*) 'Calculation complete!'
      write(7,*) 'Total calculation time: ', elapsed_time, ' seconds'	  
      
      ! ---- Output stress field ----
      call stress_fullfield_OpenMP(e3d)    

      
      ! DEALLOCATE
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
      use elas3d_mod
      implicit none

      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d
      real(dp), intent(in) :: phasemod(2,2)

      ! Locals
      integer :: i, j, k, l, m, n
      integer :: ijk, mm, nn, ii, jj, kk, ll
      integer :: nxy, i3, i8, j3, k3, n8, m3, m8
      integer :: is(nnode_fe)
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
      real(dp) :: ck(6,6), cmu(6,6), g(ngauss,ngauss,ngauss)
      real(dp) :: es(6,nnode_fe,ndof), delta(nnode_fe,ndof)
      real(dp) :: left(6), right(6)
      real(dp) :: sumval, C
      real(dp) :: x, y, z
      real(dp) :: exx, eyy, ezz, exz, eyz, exy
      real(dp), allocatable :: dk(:,:,:,:,:), b(:,:)


      allocate(dk(nphmax,nnode_fe,ndof,nnode_fe,ndof), b(ns,ndof))

      g = 0.d0;   es = 0.d0;   delta = 0.d0
      dk = 0.d0;  ck = 0.d0;   cmu = 0.d0
      left = 0.d0; right = 0.d0; b = 0.d0; C = 0.d0
      dndx = 0.d0; dndy = 0.d0; dndz = 0.d0
      e3d%cmod = 0.d0;   e3d%dk = 0.d0;   e3d%b = 0.d0
      e3d%C = 0.d0

      ! Node mapping for corners of each element
      is = (/27, 3, 2, 1, 26, 19, 18, 17/)
      nxy = nx * ny


      ! Setup bulk (ck) and shear (cmu) matrices
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

      ! Compose cmod for each phase
      do k = 1,2; do j = 1,6; do i = 1,6;
         e3d%cmod(k,i,j) = phasemod(k,1)*ck(i,j) + phasemod(k,2)*cmu(i,j)
      end do; end do; end do    

      ! Set up Simpson's integration rule weight vector
      do k3 = 1, ngauss; do j3 = 1, ngauss; do i3 = 1, ngauss;
        n = 0
        if(i3==2) n = n+1
        if(j3==2) n = n+1
        if(k3==2) n = n+1
        g(i3,j3,k3) = 4.0 ** n
      end do; end do; end do

      !------ Main Stiffness Assembly ------
      !$omp parallel do collapse(1) default(shared) private(ijk,k3,j3,i3,x,y,z,dndx,dndy,dndz,es,n8,mm,nn,ii,jj,left,right,sumval)
      do ijk = 1, nphase; do k3 = 1, ngauss; do j3 = 1, ngauss; do i3 = 1, ngauss;
        x = dble(i3-1)/2.d0
        y = dble(j3-1)/2.d0
        z = dble(k3-1)/2.d0

        ! Shape function derivatives
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

        ! Strain matrix
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

        do mm = 1, ndof
          do nn = 1, ndof
            do ii = 1, nnode_fe
              left = es(:, ii, mm)
              do jj = 1, nnode_fe
                right = es(:, jj, nn)
                sumval = dot_product(left, matmul(e3d%cmod(ijk,:,:), right))
                !$omp atomic
                dk(ijk,ii,mm,jj,nn) = dk(ijk,ii,mm,jj,nn) + g(i3,j3,k3)*sumval/216.d0
              end do
            end do
          end do
        end do

      end do; end do; end do; end do;
      !$omp end parallel do

      e3d%dk = dk ! assign back

      ! --- Set up constants for energy: b and C ---
      b = e3d%b
      C = e3d%C
	  
	  ! ---- Boundary and periodic constants ----
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy

      ! Setup delta arrays and update b, C for three directions (faces, edges, corners) 
      ! -- Begin boundary and periodic terms section --
	  
      ! -------- x = nx face --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==2 .or. i8==3 .or. i8==6 .or. i8==7) then
          delta(i8,1) = exx*nx
          delta(i8,2) = exy*nx
          delta(i8,3) = exz*nx
        end if
      end do;

      !$omp parallel do private(j,k,m,nn,mm,sumval,m3,m8)
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
              !$omp atomic
              b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- y = ny face --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==3 .or. i8==4 .or. i8==7 .or. i8==8) then
          delta(i8,1) = exy*ny
          delta(i8,2) = eyy*ny
          delta(i8,3) = eyz*ny
        end if
      end do;

      !$omp parallel do private(i,k,m,nn,mm,sumval,m3,m8)
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
              !$omp atomic
              b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- z = nz face --------
      delta = 0.d0
      do i8 = 1, nnode_fe;
        if (i8==5 .or. i8==6 .or. i8==7 .or. i8==8) then
          delta(i8,1) = exz*nz
          delta(i8,2) = eyz*nz
          delta(i8,3) = ezz*nz
        end if
      end do;

      !$omp parallel do private(i,j,m,nn,mm,sumval,m3,m8)
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
              !$omp atomic
              b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
            end do
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- x=nx, y=ny edge --------
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

      !$omp parallel do private(k,m,nn,mm,sumval,m3,m8)
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
            !$omp atomic
            b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- x=nx, z=nz edge --------
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

      !$omp parallel do private(j,m,nn,mm,sumval,m3,m8)
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
            !$omp atomic
            b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- y=ny, z=nz edge --------
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

      !$omp parallel do private(i,m,nn,mm,sumval,m3,m8)
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
            !$omp atomic
            b(e3d%ib(m,is(mm)),nn) = b(e3d%ib(m,is(mm)),nn) + sumval
          end do
        end do
      end do
      !$omp end parallel do

      ! -------- x=nx, y=ny, z=nz corner --------
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
	  ! -- End boundary and periodic terms section --

      ! Save back
      e3d%b = b
      e3d%C = C
   

      deallocate(dk, b)

      end subroutine femat_OpenMP    
!-------------------------------------------------------------------------
      subroutine energy_OpenMP(e3d, utot)
      use elas3d_mod
      implicit none

      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d
      real(dp), intent(out) :: utot

      ! Locals
      integer :: m, m3, j, n
      real(dp) :: C
      real(dp), allocatable :: gb(:,:)
      real(dp) :: gbsum

      ! Allocate locals
      allocate( gb(ns,ndof) )

      ! Unpack from struct
      C   = e3d%C

      ! Matrix operation: gb = A * u, using all local stiffness matrices
      gb(:,:) = 0.d0
      do j = 1, ndof
        !$omp parallel do private(m, n, gbsum)
        do m = 1, ns
          gbsum = 0.0d0
          do n = 1, ndof
            gbsum = gbsum                                   &
            + e3d%u(e3d%ib(m,1),n )*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,4,n) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,3,n) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,8,n) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,7,n) )   &
            + e3d%u(e3d%ib(m,2),n)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,3,n) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,7,n) )                                                                                   &
            + e3d%u(e3d%ib(m,3),n)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,2,n) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,3,n) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,7,n) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,6,n) )    &
            + e3d%u(e3d%ib(m,4),n)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,2,n) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,6,n) )                                                                                    &
            + e3d%u(e3d%ib(m,5),n)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,2,n) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,1,n) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,6,n) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,5,n) )     &
            + e3d%u(e3d%ib(m,6),n)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,1,n) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,5,n) )                                                                                    &
            + e3d%u(e3d%ib(m,7),n)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,4,n) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,1,n) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,8,n) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,5,n) )     &
            + e3d%u(e3d%ib(m,8),n)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,4,n) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,8,n) )                                                                                    &
            + e3d%u(e3d%ib(m,9),n)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,4,n) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,3,n) )                                                                                   &
            + e3d%u(e3d%ib(m,10),n)*( e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,3,n) )                                                                                                                          &
            + e3d%u(e3d%ib(m,11),n)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,3,n) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,2,n) )                                                                                  &
            + e3d%u(e3d%ib(m,12),n)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,2,n) )                                                                                                                          &
            + e3d%u(e3d%ib(m,13),n)*( e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,1,n) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,2,n) )                                                                                  &
            + e3d%u(e3d%ib(m,14),n)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,1,n) )                                                                                                                          &
            + e3d%u(e3d%ib(m,15),n)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,4,n) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,1,n) )                                                                                  &
            + e3d%u(e3d%ib(m,16),n)*( e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,4,n) )                                                                                                                          &
            + e3d%u(e3d%ib(m,17),n)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,8,n) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,7,n) )                                                                                   &
            + e3d%u(e3d%ib(m,18),n)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,7,n) )                                                                                                                          &
            + e3d%u(e3d%ib(m,19),n)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,6,n) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,7,n) )                                                                                   &
            + e3d%u(e3d%ib(m,20),n)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,6,n) )                                                                                                                           &
            + e3d%u(e3d%ib(m,21),n)*( e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,5,n) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,6,n) )                                                                                    &
            + e3d%u(e3d%ib(m,22),n)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,5,n) )                                                                                                                           &
            + e3d%u(e3d%ib(m,23),n)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,8,n) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,5,n) )                                                                                    &
            + e3d%u(e3d%ib(m,24),n)*( e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,8,n) )                                                                                                                           &
            + e3d%u(e3d%ib(m,25),n)*( e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,3,n) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,4,n) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,2,n) + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,1,n) )  &
            + e3d%u(e3d%ib(m,26),n)*( e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,7,n) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,8,n) + e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,5,n) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,6,n) )     &
            + e3d%u(e3d%ib(m,27),n)*( e3d%dk(e3d%pix(e3d%ib(m,27)),1,j,1,n) + e3d%dk(e3d%pix(e3d%ib(m,7)),2,j,2,n) + e3d%dk(e3d%pix(e3d%ib(m,6)),3,j,3,n) + e3d%dk(e3d%pix(e3d%ib(m,5)),4,j,4,n)       &
                                    + e3d%dk(e3d%pix(e3d%ib(m,25)),5,j,5,n) + e3d%dk(e3d%pix(e3d%ib(m,15)),6,j,6,n) + e3d%dk(e3d%pix(e3d%ib(m,14)),7,j,7,n) + e3d%dk(e3d%pix(e3d%ib(m,13)),8,j,8,n) )
          end do
          gb(m, j) = gbsum
        end do
        !$omp end parallel do
      end do

      ! Energy sum and gb update
      utot = C
      do m3 = 1, ndof
        !$omp parallel do reduction(+:utot)
        do m = 1, ns
          utot = utot + 0.5d0 * e3d%u(m, m3) * gb(m, m3) + e3d%b(m, m3) * e3d%u(m, m3)
          gb(m, m3) = gb(m, m3) + e3d%b(m, m3)
        end do
        !$omp end parallel do
      end do

      ! Assign back to struct
      e3d%gb = gb

      ! Deallocate locals
      deallocate( gb )

      end subroutine energy_OpenMP
!-------------------------------------------------------------------------  
      subroutine dembx_OpenMP(e3d, Lstep, gg, gg0, kkk)
      use elas3d_mod
      implicit none

      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d
      integer, intent(out) :: Lstep
      real(dp),intent(out) :: gg
      real(dp), intent(in) :: gg0
      integer, intent(in) :: kkk

      ! Locals
      integer :: ns_local, ndof_local, L1, L2
      integer :: m, ijk, i, j, n
      real(dp) :: lambda, gamma, hAh, gglast, ahsum
	  real(dp) :: rel_res
      real(dp), allocatable :: u(:,:), gb(:,:), h(:,:)
      real(dp), allocatable :: Ah(:,:)
      integer :: ibl(27), pl(27)
      
      ns_local   = ns
      ndof_local = ndof      

      ! Allocate locals
      allocate( u(ns,ndof), gb(ns,ndof), h(ns,ndof) )
      allocate( Ah(ns,ndof) )

      ! Unpack from struct
      u   = e3d%u
      gb  = e3d%gb
      h   = e3d%h


      ! Initialize conjugate direction vector
      if (kkk == 1) h(:,:) = gb(:,:)

      Lstep = 0
      do ijk = 1, ldemb
        Lstep = Lstep + 1

        ! Zero out Ah
        Ah = 0.d0

        !--- Blocking over ns for better cache locality ---
        do L1 = 1, ns_local, block_size
          L2 = min(L1 + block_size - 1, ns_local)        
          
          ! Matrix multiply: Ah = A * h
		  !$omp parallel do private(m, n, j, ibl, pl, ahsum) schedule(guided)
          do m = L1, L2
            ! Caching
            ibl = e3d%ib(m,1:27)     ! All 27 neighbors
            pl  = e3d%pix(ibl(1:27)) ! Phase index for all neighbors
            do j = 1, ndof_local
              ! Calculation
              ahsum = 0.d0
              !$omp simd reduction(+:ahsum)
              do n = 1, ndof_local
                ahsum = ahsum &
               + h(ibl(1),n )*( e3d%dk(pl(27),1,j,4,n) + e3d%dk(pl(7),2,j,3,n) + e3d%dk(pl(25),5,j,8,n) + e3d%dk(pl(15),6,j,7,n) ) &
               + h(ibl(2),n)*( e3d%dk(pl(27),1,j,3,n) + e3d%dk(pl(25),5,j,7,n) ) &
               + h(ibl(3),n)*( e3d%dk(pl(27),1,j,2,n) + e3d%dk(pl(5),4,j,3,n) + e3d%dk(pl(13),8,j,7,n) + e3d%dk(pl(25),5,j,6,n) ) &
               + h(ibl(4),n)*( e3d%dk(pl(5),4,j,2,n) + e3d%dk(pl(13),8,j,6,n) ) &
               + h(ibl(5),n)*( e3d%dk(pl(6),3,j,2,n) + e3d%dk(pl(5),4,j,1,n) + e3d%dk(pl(14),7,j,6,n) + e3d%dk(pl(13),8,j,5,n) ) &
               + h(ibl(6),n)*( e3d%dk(pl(6),3,j,1,n) + e3d%dk(pl(14),7,j,5,n) ) &
               + h(ibl(7),n)*( e3d%dk(pl(6),3,j,4,n) + e3d%dk(pl(7),2,j,1,n) + e3d%dk(pl(14),7,j,8,n) + e3d%dk(pl(15),6,j,5,n) ) &
               + h(ibl(8),n)*( e3d%dk(pl(7),2,j,4,n) + e3d%dk(pl(15),6,j,8,n) ) &
               + h(ibl(9),n)*( e3d%dk(pl(25),5,j,4,n) + e3d%dk(pl(15),6,j,3,n) ) &
               + h(ibl(10),n)*( e3d%dk(pl(25),5,j,3,n) ) &
               + h(ibl(11),n)*( e3d%dk(pl(13),8,j,3,n) + e3d%dk(pl(25),5,j,2,n) ) &
               + h(ibl(12),n)*( e3d%dk(pl(13),8,j,2,n) ) &
               + h(ibl(13),n)*( e3d%dk(pl(13),8,j,1,n) + e3d%dk(pl(14),7,j,2,n) ) &
               + h(ibl(14),n)*( e3d%dk(pl(14),7,j,1,n) ) &
               + h(ibl(15),n)*( e3d%dk(pl(14),7,j,4,n) + e3d%dk(pl(15),6,j,1,n) ) &
               + h(ibl(16),n)*( e3d%dk(pl(15),6,j,4,n) ) &
               + h(ibl(17),n)*( e3d%dk(pl(27),1,j,8,n) + e3d%dk(pl(7),2,j,7,n) ) &
               + h(ibl(18),n)*( e3d%dk(pl(27),1,j,7,n) ) &
               + h(ibl(19),n)*( e3d%dk(pl(27),1,j,6,n) + e3d%dk(pl(5),4,j,7,n) ) &
               + h(ibl(20),n)*( e3d%dk(pl(5),4,j,6,n) ) &
               + h(ibl(21),n)*( e3d%dk(pl(5),4,j,5,n) + e3d%dk(pl(6),3,j,6,n) ) &
               + h(ibl(22),n)*( e3d%dk(pl(6),3,j,5,n) ) &
               + h(ibl(23),n)*( e3d%dk(pl(6),3,j,8,n) + e3d%dk(pl(7),2,j,5,n) ) &
               + h(ibl(24),n)*( e3d%dk(pl(7),2,j,8,n) ) &
               + h(ibl(25),n)*( e3d%dk(pl(14),7,j,3,n) + e3d%dk(pl(13),8,j,4,n) + e3d%dk(pl(15),6,j,2,n) + e3d%dk(pl(25),5,j,1,n) ) &
               + h(ibl(26),n)*( e3d%dk(pl(6),3,j,7,n) + e3d%dk(pl(5),4,j,8,n) + e3d%dk(pl(27),1,j,5,n) + e3d%dk(pl(7),2,j,6,n) ) &
               + h(ibl(27),n)*( e3d%dk(pl(27),1,j,1,n) + e3d%dk(pl(7),2,j,2,n) + e3d%dk(pl(6),3,j,3,n) + e3d%dk(pl(5),4,j,4,n) &
                                + e3d%dk(pl(25),5,j,5,n) + e3d%dk(pl(15),6,j,6,n) + e3d%dk(pl(14),7,j,7,n) + e3d%dk(pl(13),8,j,8,n) )
              end do
              Ah(m, j) = ahsum
            end do
          end do
          !$omp end parallel do
        end do ! blocks          

        ! Compute hAh = sum(h .* Ah)
        ! Parallel reduction for hAh and
		! Standard conjugate gradient update(gg)
        hAh = 0.d0
		gg = 0.d0
        !$omp parallel do reduction(+:hAh, gg) collapse(2)
        do i = 1, ns_local; do j = 1, ndof_local
          hAh = hAh + h(i,j) * Ah(i,j)
		  gg = gg + gb(i,j)*gb(i,j)
        end do; end do;
        !$omp end parallel do        


        if(dabs(hAh) < 1d-14) then
          print *,"WARNING: hAh is nearly zero, breaking CG step."
          exit
        end if
    
		! Compute Step Size        
        lambda = gg / hAh

        ! Solution/Residual update
        !$omp parallel do collapse(2)
        do i = 1, ns; do j = 1, ndof
          u(i,j) = u(i,j) - lambda * h(i,j)
          gb(i,j) = gb(i,j) - lambda * Ah(i,j)
        end do; end do;
        !$omp end parallel do


        ! Update new norm
        gglast = gg
        gg = 0.d0
        !$omp parallel do reduction(+:gg) collapse(2)
        do i = 1, ns; do j = 1, ndof
          gg = gg + gb(i,j) * gb(i,j)
        end do; end do;
        !$omp end parallel do

        ! Exit early if converged
        rel_res = sqrt(gg/gg0)
        if (rel_res < tol) exit

        ! Conjugate direction update (Fletcher–Reeves)
        gamma = gg / gglast
		!$omp parallel do collapse(2)
        do i = 1, ns; do j = 1, ndof
          h(i,j) = gb(i,j) + gamma * h(i,j)
        end do; end do;
        !$omp end parallel do
        
        
        print *, 'Number of conjugate steps_local iteration =', Lstep
        print *, 'rel. residual norm=', rel_res      
        
      end do

      ! Save back to struct fields
      e3d%u  = u
      e3d%gb = gb
      e3d%h  = h


      deallocate( u, gb, h, Ah )

      end subroutine dembx_OpenMP
!-------------------------------------------------------------------------
      subroutine stress_OpenMP(e3d)
      use elas3d_mod
      implicit none

      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d

      ! Locals
      integer :: i, j, k, m, nxy, mm, n8, n3, n
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
      real(dp) :: uu(nnode_fe,ndof), es(6,nnode_fe,ndof)
      real(dp) :: str11, str22, str33, str13, str23, str12
      real(dp) :: s11, s22, s33, s13, s23, s12
      real(dp) :: exx, eyy, ezz, exz, eyz, exy
      ! Output fields	  
      real(dp) :: strxx, stryy, strzz, strxz, stryz, strxy
      real(dp) :: sxx, syy, szz, sxz, syz, sxy
      

      ! For local index
      integer :: ibl(27)
      integer :: m_local, k_local, j_local, i_local, pixm
      real(dp) :: str11_t, str22_t, str33_t, str13_t, str23_t, str12_t
      real(dp) :: s11_t, s22_t, s33_t, s13_t, s23_t, s12_t


      ! Struct unpack
      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy

      nxy = nx * ny

      ! --- Reference shape-function derivatives ---
      dndx = [-0.25d0, 0.25d0, 0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0, -0.25d0]
      dndy = [-0.25d0, -0.25d0, 0.25d0, 0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0]
      dndz = [-0.25d0, -0.25d0, -0.25d0, -0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0]

      ! Build strain matrix
      es = 0.d0
      do n = 1, nnode_fe
        es(1,n,1) = dndx(n)
        es(2,n,2) = dndy(n)
        es(3,n,3) = dndz(n)
        es(4,n,1) = dndz(n)
        es(4,n,3) = dndx(n)
        es(5,n,2) = dndz(n)
        es(5,n,3) = dndy(n)
        es(6,n,1) = dndy(n)
        es(6,n,2) = dndx(n)
      end do


      strxx = 0.d0; stryy = 0.d0; strzz = 0.d0
      strxz = 0.d0; stryz = 0.d0; strxy = 0.d0
      sxx = 0.d0;   syy = 0.d0;   szz = 0.d0
      sxz = 0.d0;   syz = 0.d0;   sxy = 0.d0

      !$omp parallel default(shared) private(k_local,j_local,i_local,m_local,mm,n8,n3,n,ibl,uu,str11_t,str22_t,str33_t,str13_t,str23_t,str12_t, s11_t,s22_t,s33_t,s13_t,s23_t,s12_t, pixm) reduction(+:strxx, stryy, strzz, strxz, stryz, strxy, sxx, syy, szz, sxz, syz, sxy)
      do k_local = 1, nz
        do j_local = 1, ny
          do i_local = 1, nx
            m_local = (k_local-1)*nxy + (j_local-1)*nx + i_local

            ! Cache ib for m_local
            ibl( 1) = e3d%ib(m_local,  1); ibl( 2) = e3d%ib(m_local,  2); 
            ibl( 3) = e3d%ib(m_local,  3); 
            ibl(17) = e3d%ib(m_local, 17); ibl(18) = e3d%ib(m_local, 18); 
            ibl(19) = e3d%ib(m_local, 19); ibl(26) = e3d%ib(m_local, 26);

            ! Load elements for 8-vector using periodic boundary conditions
            uu = 0.d0
            do mm = 1, ndof
              uu(1,mm) = e3d%u(m_local,mm)
              uu(2,mm) = e3d%u(ibl(3),mm)
              uu(3,mm) = e3d%u(ibl(2),mm)
              uu(4,mm) = e3d%u(ibl(1),mm)
              uu(5,mm) = e3d%u(ibl(26),mm)
              uu(6,mm) = e3d%u(ibl(19),mm)
              uu(7,mm) = e3d%u(ibl(18),mm)
              uu(8,mm) = e3d%u(ibl(17),mm)
            end do

            ! Periodic boundary corrections
            if(i_local == nx) then
              uu(2,1) = uu(2,1) + exx*nx
              uu(3,1) = uu(3,1) + exx*nx
              uu(6,1) = uu(6,1) + exx*nx
              uu(7,1) = uu(7,1) + exx*nx
              uu(2,2) = uu(2,2) + exy*nx
              uu(3,2) = uu(3,2) + exy*nx
              uu(6,2) = uu(6,2) + exy*nx
              uu(7,2) = uu(7,2) + exy*nx
              uu(2,3) = uu(2,3) + exz*nx
              uu(3,3) = uu(3,3) + exz*nx
              uu(6,3) = uu(6,3) + exz*nx
              uu(7,3) = uu(7,3) + exz*nx
            end if
            if(j_local == ny) then
              uu(3,1) = uu(3,1) + exy*ny
              uu(4,1) = uu(4,1) + exy*ny
              uu(7,1) = uu(7,1) + exy*ny
              uu(8,1) = uu(8,1) + exy*ny
              uu(3,2) = uu(3,2) + eyy*ny
              uu(4,2) = uu(4,2) + eyy*ny
              uu(7,2) = uu(7,2) + eyy*ny
              uu(8,2) = uu(8,2) + eyy*ny
              uu(3,3) = uu(3,3) + eyz*ny
              uu(4,3) = uu(4,3) + eyz*ny
              uu(7,3) = uu(7,3) + eyz*ny
              uu(8,3) = uu(8,3) + eyz*ny
            end if
            if(k_local == nz) then
              uu(5,1) = uu(5,1) + exz*nz
              uu(6,1) = uu(6,1) + exz*nz
              uu(7,1) = uu(7,1) + exz*nz
              uu(8,1) = uu(8,1) + exz*nz
              uu(5,2) = uu(5,2) + eyz*nz
              uu(6,2) = uu(6,2) + eyz*nz
              uu(7,2) = uu(7,2) + eyz*nz
              uu(8,2) = uu(8,2) + eyz*nz
              uu(5,3) = uu(5,3) + ezz*nz
              uu(6,3) = uu(6,3) + ezz*nz
              uu(7,3) = uu(7,3) + ezz*nz
              uu(8,3) = uu(8,3) + ezz*nz
            end if

            ! Local accumulators
            str11_t = 0.d0; str22_t = 0.d0; str33_t = 0.d0
            str13_t = 0.d0; str23_t = 0.d0; str12_t = 0.d0
            s11_t = 0.d0;   s22_t = 0.d0;   s33_t = 0.d0
            s13_t = 0.d0;   s23_t = 0.d0;   s12_t = 0.d0

            pixm = e3d%pix(m_local)
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
      !$omp end parallel

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
      use elas3d_mod
      implicit none
      
      ! Arguments
      type(elas3d_data_type), intent(in) :: e3d
      real(dp), intent(out) :: prob(nphase)
      
      ! Locals
      integer :: i, m

      ! Zero the phase counters
      prob(:) = 0.d0

      ! Count the occurrences of each phase
      do m = 1, ns
        if(e3d%pix(m) >= 1 .and. e3d%pix(m) <= nphase) then
          prob(e3d%pix(m)) = prob(e3d%pix(m)) + 1.d0
        end if
      end do

      ! Convert counts to volume fractions
      prob(:) = prob(:) / dble(ns) 

      end subroutine assig
!-------------------------------------------------------------------------
      subroutine stress_fullfield_OpenMP(e3d) 
      use elas3d_mod
      use hdf5
      implicit none
  
      ! Arguments
      type(elas3d_data_type), intent(inout) :: e3d

      ! Locals
      integer :: i, j, k, m, nxy, mm, n8, n3, n
      real(dp) :: dndx(nnode_fe), dndy(nnode_fe), dndz(nnode_fe)
      real(dp) :: uu(nnode_fe,ndof), es(6,nnode_fe,ndof)
      real(dp) :: eps(6)
      real(dp) :: exx, eyy, ezz, exz, eyz, exy      
      
      ! For local index
      integer :: ibl(27)
      integer :: m_l, k_l, j_l, i_l      

      ! HDF5 variables
      integer(HID_T)    :: file_id, dset_id, space_id
      integer(HSIZE_T)  :: dims(1)
      integer           :: hdferr   


      exx = e3d%exx; eyy = e3d%eyy; ezz = e3d%ezz
      exz = e3d%exz; eyz = e3d%eyz; exy = e3d%exy  
  
      ! Physical sizes
      nxy = nx * ny  

      ! --- Brick basis functions
      dndx = [-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0,-0.25d0]
      dndy = [-0.25d0,-0.25d0, 0.25d0, 0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0]
      dndz = [-0.25d0,-0.25d0,-0.25d0,-0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0]
      
      ! Assemble engineering strain Voigt matrix	  
      es = 0.d0
      do n = 1,nnode_fe
        es(1,n,1) = dndx(n)
        es(2,n,2) = dndy(n)
        es(3,n,3) = dndz(n)
        es(4,n,1) = dndz(n);    es(4,n,3) = dndx(n)
        es(5,n,2) = dndz(n);    es(5,n,3) = dndy(n)
        es(6,n,1) = dndy(n);    es(6,n,2) = dndx(n)
      end do

      
      ! Main per-voxel calculation
      !$omp parallel do default(shared) private(k_l,j_l,i_l,m_l,mm,n8,n3,n,ibl,uu,eps)
      do k_l = 1, nz; do j_l = 1, ny; do i_l = 1, nx;
        m_l = (k_l-1)*nxy + (j_l-1)*nx + i_l

        ! Cache ib for m_l
        ibl( 1) = e3d%ib(m_l,  1); 
		ibl( 2) = e3d%ib(m_l,  2); 
		ibl( 3) = e3d%ib(m_l,  3); 
        ibl(17) = e3d%ib(m_l, 17); 
		ibl(18) = e3d%ib(m_l, 18); 
		ibl(19) = e3d%ib(m_l, 19); 
		ibl(26) = e3d%ib(m_l, 26);
        
        ! Local displacements
        uu = 0.d0
        do mm = 1, ndof
          uu(1,mm) = e3d%u(m_l,mm)
          uu(2,mm) = e3d%u(ibl(3),mm)
          uu(3,mm) = e3d%u(ibl(2),mm)
          uu(4,mm) = e3d%u(ibl(1),mm)
          uu(5,mm) = e3d%u(ibl(26),mm)
          uu(6,mm) = e3d%u(ibl(19),mm)
          uu(7,mm) = e3d%u(ibl(18),mm)
          uu(8,mm) = e3d%u(ibl(17),mm)
        end do
    
        ! Periodic boundary conditions
        if (i_l == nx) then
          uu(2,1) = uu(2,1) + exx*nx
          uu(3,1) = uu(3,1) + exx*nx
          uu(6,1) = uu(6,1) + exx*nx
          uu(7,1) = uu(7,1) + exx*nx
          uu(2,2) = uu(2,2) + exy*nx
          uu(3,2) = uu(3,2) + exy*nx
          uu(6,2) = uu(6,2) + exy*nx
          uu(7,2) = uu(7,2) + exy*nx
          uu(2,3) = uu(2,3) + exz*nx
          uu(3,3) = uu(3,3) + exz*nx
          uu(6,3) = uu(6,3) + exz*nx
          uu(7,3) = uu(7,3) + exz*nx
        end if
        if (j_l == ny) then
          uu(3,1) = uu(3,1) + exy*ny
          uu(4,1) = uu(4,1) + exy*ny
          uu(7,1) = uu(7,1) + exy*ny
          uu(8,1) = uu(8,1) + exy*ny
          uu(3,2) = uu(3,2) + eyy*ny
          uu(4,2) = uu(4,2) + eyy*ny
          uu(7,2) = uu(7,2) + eyy*ny
          uu(8,2) = uu(8,2) + eyy*ny
          uu(3,3) = uu(3,3) + eyz*ny
          uu(4,3) = uu(4,3) + eyz*ny
          uu(7,3) = uu(7,3) + eyz*ny
          uu(8,3) = uu(8,3) + eyz*ny
        end if
        if (k_l == nz) then
          uu(5,1) = uu(5,1) + exz*nz
          uu(6,1) = uu(6,1) + exz*nz
          uu(7,1) = uu(7,1) + exz*nz
          uu(8,1) = uu(8,1) + exz*nz
          uu(5,2) = uu(5,2) + eyz*nz
          uu(6,2) = uu(6,2) + eyz*nz
          uu(7,2) = uu(7,2) + eyz*nz
          uu(8,2) = uu(8,2) + eyz*nz
          uu(5,3) = uu(5,3) + ezz*nz
          uu(6,3) = uu(6,3) + ezz*nz
          uu(7,3) = uu(7,3) + ezz*nz
          uu(8,3) = uu(8,3) + ezz*nz
        end if         

        
        !--- Compute engineering strains
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

        ! Cauchy stress
        e3d%stress_field(m_l,:) = matmul(e3d%cmod(e3d%pix(m_l),1:6,1:6), eps)

        
      end do; end do; end do;
      !$omp end parallel do

      
      ! Von Mises stress
      !$omp parallel do
      do m = 1, ns
        e3d%vm(m) = sqrt(0.5d0*((e3d%stress_field(m,1)-e3d%stress_field(m,2))**2 + &
                                 (e3d%stress_field(m,2)-e3d%stress_field(m,3))**2 + &
                                 (e3d%stress_field(m,3)-e3d%stress_field(m,1))**2 + &
                                  6.0d0*(e3d%stress_field(m,4)**2 + &
                                         e3d%stress_field(m,5)**2 + &
										 e3d%stress_field(m,6)**2)))
      end do
      !$omp end parallel do


      

      ! ========== HDF5 EXPORT ==========
      dims(1) = ns
      call h5open_f(hdferr)
      ! Create file
      call h5fcreate_f('fullfield_poly.h5', H5F_ACC_TRUNC_F, file_id, hdferr)
      ! Create dataspace
      call h5screate_simple_f(1, dims, space_id, hdferr)
	  ! Create dataset (integers), Save pix(m)
      call h5dcreate_f(file_id, "pix", H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
	  call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, e3d%pix, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
	  ! Create dataset (doubles), Save vm
      call h5dcreate_f(file_id, "vm", H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, e3d%vm, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
      
      call h5sclose_f(space_id, hdferr)
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

      
      end subroutine stress_fullfield_OpenMP    
!-------------------------------------------------------------------------    
      subroutine ppixel_hdf5(e3d)
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


      ! === Read pix (microstructure phase assignment) ===
      call h5open_f(hdferr)
      call h5fopen_f('input_structure_poly.h5', H5F_ACC_RDONLY_F, file_id, hdferr)
	  if (hdferr /= 0) stop 'Error opening HDF5 file'	 	  
      call h5dopen_f(file_id, 'pix', dset_id, hdferr)
      call h5dget_space_f(dset_id, space_id, hdferr)

      
      ! Get dimensions to verify size (rank)
      call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
      allocate(dims(rank), maxdims(rank))
      call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)

	  total_elems = product(dims)	  

      if (total_elems /= ns) then
        print *, 'ERROR: HDF5 pix ns not matching F90 ns!',  product(dims), ns
        stop
      end if

      ! Read Data (Use H5S_ALL_F to read the whole dataset)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, e3d%pix, dims, hdferr)
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      deallocate(dims, maxdims)

      ! === Read orientation field ===
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
	  
	  ! Read Data	  
      odims(1) = size(e3d%orientation)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, e3d%orientation, odims, hdferr)
	  
	  ! Cleanup 'orientation' handles		  
      call h5dclose_f(dset_id, hdferr)
      call h5sclose_f(space_id, hdferr)
      deallocate(dims, maxdims)


      ! Finished
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

      ! Checks
      if (minval(e3d%pix) < 1 .or. maxval(e3d%pix) > nphase) then
        print *, "ERROR: Phase label in pix out of bounds"
        stop
      end if

      end subroutine ppixel_hdf5      