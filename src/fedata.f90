module fedata
    use types
    implicit none    
    save

    !! This module is used to define global values and matrices
    
    integer :: ne
        !! Total number of elements
    integer :: nn
        !! Total number of nodes
    integer :: nb
        !! Number of boundary conditions
    integer :: np
        !! Number of loads
    integer :: nm
        !! Number of materials
    integer :: neqn
        !! Number of equations
    integer :: bw
        !! Bandwidth of system matrix
    
    real(wp), dimension(:,:), allocatable :: x
        !! Nodal coordinates:
        !! * _i_-th row: coordinates for node _i_
        !! * column 1: \(x\)-coordinate
        !! * column 2: \(y\)-coordinate
        
    type element_def
        !! The "element" data type, which is used to
        !! describe a single element of a structure
        integer, dimension(4) :: ix
            !! Node list
        integer :: id
            !! Element id:
            !! * 1: truss element
            !! * 2: continuum element
        integer :: mat
            !! Identification number for material used by this element
        integer :: numnode
            !! Number of nodes on this element
    end type element_def
    type(element_def), dimension(:), allocatable :: element
        !! The elements of a structure

    type matprop
        !! The "material" data type, which is used to
        !! describe a material. Note that a given material
        !! might not use all of the properties defined below
        real(wp) :: young
            !! Young's Modulus
        real(wp) :: nu
            !! Poisson's Ratio
        real(wp) :: thk
            !! The thickness of the structure where this material is used
        real(wp) :: area
            !! The cross-sectional area of the structure where this material is used
        real(wp) :: dens
            !! Density
        real(wp) :: youngy
            !! Young's Modulus in the transverse direction (used for orthotropic materials)
        real(wp) :: shear
            !! Shear Modulus (used for orthotropic materials)
    end type matprop
    type(matprop), dimension(:), allocatable :: mprop
        !! The materials of a structure
        
    ! Boundary conditions
    real(wp), dimension(:,:), allocatable :: bound
        !! Boundary conditions
        !!
        !! * _i_-th row: boundary condition number _i_
        !! * column 1: index of node where boundary condition is applied
        !! * column 2: degree of freedom _i_ affected by boundary condition
        !!     * _i_ = 1: \(x\)-component
        !!     * _i_ = 2: \(y\)-component
        !! * column 3: prescribed value for the displacement
    real(wp), dimension(:,:), allocatable :: loads
        !! External loading
        !!
        !! * _i_-th row: external load number _i_
        !! * column 1: type _j_ of external loading
        !!     * _j_ = 1: point load
        !!     * _j_ = 2: pressure load
        !!
        !! The meaning of columns 2 - 4 depend on the type of external loading:
        !!
        !! * For point loads:
        !!     * column 2: index of node where load is applied
        !!     * column 3: degree of freedom _k_ affected by load
        !!         - _k_ = 1: \(x\)-component
        !!         - _k_ = 2: \(y\)-component
        !!      * column 4: value of applied force
        !!
        !! * For surface loads:
        !!     * column 2: index of element where pressure load is applied
        !!     * column 3: index of element face where pressure is applied
        !!     * column 4: value of applied pressure (always normal to surface)
    real(wp) :: accel(2)
        !! Acceleration forces
        !!
        !! * index 1: acceleration along \(x\)-component
        !! * index 2: acceleration along \(y\)-component
        
    ! Working arrays:
    real(wp), dimension(:,:), allocatable :: kmatrix
        !! Stiffness matrix
    real(wp), dimension(:,:), allocatable :: strain
        !! Strains at different places in the structure 
        !!
        !! * _i_-th row: strain in element _i_  
        !! * column 1: \(\epsilon_{11}\)  
        !! * column 2: \(\epsilon_{22}\)  
        !! * column 3: \(\epsilon_{12}\)  
    real(wp), dimension(:,:), allocatable :: stress
        !! Stresses at different places in the structure 
        !!
        !! * _i_-th row: stress in element _i_  
        !! * column 1: \(\sigma_{11}\)  
        !! * column 2: \(\sigma_{22}\)  
        !! * column 3: \(\sigma_{12}\)      
    real(wp), dimension(:),   allocatable :: p
        !! Force vector
    real(wp), dimension(:),   allocatable :: d
        !! Displacement vector
    real(wp), dimension(:,:),   allocatable :: d_eig
    	!! Eigenvector
    real(wp), dimension(:,:),   allocatable :: d_all
    	!! History of displacements
    real(wp), dimension(:), allocatable :: vm_stress
    real(wp), dimension(:), allocatable :: psi
    	!! Von mises Stress and direction (psi)

    ! Input/Output
    character(len=50) :: filename
        !! Name of input file
!    character(len=5), parameter :: plotdevice = '/XWIN'
    character(len=3), parameter :: plotdevice = '/SS'    
!    character(len=3), parameter :: plotdevice = '/GW'
        !! Platform (or graphics package) dependent parameter for plotting
        !!
        !! * Linux               : '/XWIN'  
        !! * Windows/Plato       : '/SS'
        !! * Windows/Code::Blocks: '/GW'

    ! Analysis type
    character(len=20) :: antype
        !! The type of analysis to be carried out
        !!
        !! STATIC
        !! : Static linear analysis (default)
        !!
        !! STATIC_NL
        !! : Static geometrically nonlinear analysis (HINT: not implemented yet)
        !!
        !! MODAL
        !! : Analysis of eigenfrequencies and mode shapes (HINT: not implemented yet)
        !!
        !! ANGLE
        !! : Optimization of fiber orientation of anisotropic material (HINT: not implemented yet)
        !!
        !! TRANS
        !! : Transient analysis (HINT: not implemented yet)

    ! Parameters
    integer, parameter :: mdim = 8
        !! Maximum number of element degrees-of-freedom
    integer, parameter :: no_gauss_points_k = 2
    	!! Number of points for gauss numeric integration of stiffness
    integer, parameter :: no_gauss_points_s = 1
    	!! Number of points for gauss numeric integration of strain
    integer, parameter :: no_eigenvalues = 5
    	!! Number of eigenvalues to be found

        
    real(wp), parameter :: scale_def = 1.0_wp
        !! Scale deformations when plotting
    real(wp), parameter :: scale_vec = 1.0_wp
        !! Scale length of vectors for plotting
	logical, parameter :: plot_stress = .false.
    	!! Plot regular stress
    logical, parameter :: plot_vmstress = .true.
    	!! Plot Von Mises Stress

        
    real(wp), parameter :: scale_thk = 0.1_wp
        !! Scale thickness of lines
    real(wp), parameter :: element_evaluation_x = 0.0_wp
    real(wp), parameter :: element_evaluation_y = 0.0_wp
    	!! Set point of evaluation for continuum elements (only used for rectangular elements)
    real(wp), parameter :: estop = 10e-12_wp
		!! Set presicion for for inverse iteration in modal analysis
        
    logical, parameter :: isoparametric = .true.
    	!! Indicate whether the  element is rectangular (false) or isoparametric (true)
	logical, parameter :: banded = .false.
        !! Indicate whether the system matrix is in banded form or not (full matrix)
    logical, parameter :: penalty = .false.
        !! Indicate whether boundary conditions are imposed by the penalty method or not (zero-one method)
    real(wp), parameter :: penalty_fac = 1.0e10_wp
        !! Scaling factor for boundary conditions imposed by the penalty method

end module fedata
