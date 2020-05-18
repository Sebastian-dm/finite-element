module fea

    !! This module contains procedures that are common to FEM-analysis
    
    implicit none
    save
    private
    public :: displ, initial, buildload, buildstiff, enforce, recover, eigen, trans

contains
!
!--------------------------------------------------------------------------------------------------
!
    subroutine initial
        
        !! This subroutine is mainly used to allocate vectors and matrices
        
        use fedata
        use link1
        use plane42
        use plane42rect
               
		integer :: e, bwe
        

        ! This subroutine computes the number of global equation,
        ! half bandwidth, etc and allocates global arrays.

        ! Calculate number of equations
        neqn = 2*nn

        if (.not. banded) then
            allocate (kmatrix(neqn, neqn))
        else
          	bw = 0
          	do e = 1, ne
            	! Element bandwidth is found from element connectivity (max and min global node index) converted with idof
            	bwe = maxval(element(e)%ix)*2 - (minval(element(e)%ix)*2-1) + 1
                ! update global bandwidth value if element bandwidth is larger
                if (bwe > bw) then
                    bw = bwe
                end if
!               print *, 'Element:',e,'  bandwidth:',bwe
        	end do
            print *, 'Global bandwidth: ',bw
            allocate (kmatrix(bw, neqn))
       end if
       
       allocate (p(neqn), d(neqn), vm_stress(ne), psi(ne))
       allocate (strain(ne, 3), stress(ne, 3))

       select case(antype)
           case ('modal')
           		allocate(D_eig(neqn,no_eigenvalues))
                D_eig = 0
       end select
    
       ! Initial d stress and strain
       d = 0
       strain = 0
       stress = 0
       vm_stress(:) = 0
       psi(:) = 0        
    end subroutine initial
!
!--------------------------------------------------------------------------------------------------
!
    subroutine displ

        !! This subroutine calculates displacements

        use fedata
        use numeth
        use processor
        
        integer :: e
        real(wp), dimension(:), allocatable :: plotval, plotval_2
        character(len=16) :: plottitle

		call stopwatch('star')
        
        ! Build load-vector
        call buildload

        ! Build stiffness matrix
        call buildstiff

        ! Remove rigid body modes
        call enforce(kmatrix)

        ! Transfer results to make sure p is not overwritten
        d(1:neqn) = p(1:neqn)

        if (.not. banded) then
            ! Factor stiffness matrix
            call factor(kmatrix)
            ! Solve for displacement vector
            call solve(kmatrix, d)
        else    
            ! Factor stiffness matrix
            call bfactor(kmatrix)
            ! Solve for displacement vector, (d is overwritten and returned)
            call bsolve(kmatrix, d)
        end if

        call stopwatch('stop')
		
        
        ! Recover stress
        call recover
                
        ! Output results
        call output

        ! Plot deford shape
        call plot( undeformed + deformed )
		
		! Plot regular stress
		if (plot_stress) then
            ! Store plot values for regular stress
            allocate (plotval(ne))
            do e = 1, ne
                if (element(e)%id == 1) then
                    plotval(e) = stress(e,1)
                else if (element(e)%id == 2) then
                    plotval(e) = stress(e,1)
                end if
            end do
            plottitle = "Stress"
            call plot( elements, eval=plotval, color=.true. , title="Stress", legend=.true.)
        end if

        ! Plot Von Mises Stress
        if (plot_vmstress) then
            allocate (plotval_2(ne))
            do e = 1, ne
                if (element(e)%id == 1) then
                    plotval_2(e) = vm_stress(e)
                else if (element(e)%id == 2) then
                    plotval_2(e) = vm_stress(e)
                end if
            end do
            call plot( elements, eval=plotval_2, color=.true. , title="Von Mises Stress", legend=.true.)
        end if

		! plot to matlab
		!call plot( elements, eval=psi, color=.true. , title=title, legend=.true., device=Matlab )
                
    end subroutine displ

!
!--------------------------------------------------------------------------------------------------
!
	subroutine trans
		
		!! This subroutine calculates the transient response of a structure

        use fedata
        use numeth
        use processor
        use plane42

        character(len=500) :: outfile_path
        character(len=20) :: t

        ! Timestep variables
        integer :: no_timestep
        real(wp) :: timestep
		
		! Declare variables for building of transmatrix
        real(wp) :: me(mdim,mdim), lme(mdim), xe(mdim), cebc(mdim,mdim)
        real(wp) :: dens, thk, transmatrix(neqn,neqn), mmatrix(neqn,neqn), cmatrix(neqn,neqn)
        integer :: edof(mdim), nen, idof, n
        
        
        ! Declare rest of variables
        integer :: n, e, i, j 
        real(wp) :: MD(neqn), MD_old(neqn), KD(neqn), KD_old(neqn), CD(neqn), CD_old(neqn)
        real(wp) :: alpha, beta, timestep, force_scale
        logical :: lumped_mass, half_step, transient_force, matlab_output, boundary_c
		
		
		! User Parameters
        alpha = 0.0
        beta = 0.0
        lumped_mass = .true.
		half_step = .false.
        transient_force = .false.
        matlab_output = .false.
        boundary_c = .false.
		timestep = 0.3
		no_timestep = 3000
        force_scale = 1.0
		
		print '(A,f20.4)',"Timestep:", timestep
        print *, 'no_timestep:',no_timestep
        print *, 'alpha:',alpha
        print *, 'beta:',beta
        print *, 'Lumped:',lumped_mass
		print *, 'Half step:',half_step
        

        ! Allocate space for all historical displacements
        allocate(D_all(neqn,no_timestep))

        ! Build load-vector
        call buildload
        
		! scale load (P)
        P = P*force_scale

        !Build left side K
        call buildstiff
        !Build left side M
        mmatrix = 0
        cmatrix = 0
        do e = 1, ne
            ! Find element coordinates
            nen = element(e)%numnode
            ! Find element dofs
            do n = 1, nen
                 xe(2*n-1) = x(element(e)%ix(n),1)
                 xe(2*n  ) = x(element(e)%ix(n),2)
                 edof(2*n-1) = 2 * element(e)%ix(n) - 1  
                 edof(2*n)   = 2 * element(e)%ix(n)
            end do
            ! Load element properties
            thk = mprop(element(e)%mat)%thk ! thk - element thickness
            dens = mprop(element(e)%mat)%dens ! dens - element Density
            select case (lumped_mass)
            	case(.true.)
                	! Create element mass
        			call plane42_lme(xe, thk, dens, lme)          
            		! Add element mass to global mass
                    do i = 1,nen*2
                            mmatrix(edof(i),edof(i)) = mmatrix(edof(i),edof(i)) + lme(i)
                    end do
                case(.false.)
					! Create element mass
        			call plane42_me(xe, thk, dens, me)          
            		! Add element mass to global mass
                    do i = 1,nen*2
                        do j = 1,nen*2
                            mmatrix(edof(i),edof(j)) = mmatrix(edof(i),edof(j)) + me(i,j)
                        end do
                    end do
            end select
        end do
        !Build left side C
        select case (boundary_c)
          	case(.false.)
            	cmatrix = alpha*mmatrix + beta*kmatrix
         	case(.true.)
            	do e = 1, ne
                	! Find element dofs
                    nen = element(e)%numnode
                    do n = 1, nen
                         edof(2*n-1) = 2 * element(e)%ix(n) - 1  
                         edof(2*n)   = 2 * element(e)%ix(n)
                    end do
                    ! Build element damping matrix
                    cebc = 0
                	call plane42_cebc(e, cebc)
                    
                    ! Assemble global damping matrix
                    cmatrix(edof,edof) = cmatrix(edof,edof) + cebc
                    !print *, "element: ",e
                    !print *, "cebc"
                    !print *, cebc
                    !print*, "cmatrix(edof,edof)",cmatrix(edof,edof)
                end do
        end select
        
		!Build left side (transmatrix) from K, M and C      
		select case (half_step)
          	case(.false.)
				transmatrix = 1.0_wp/(timestep**2.0_wp)*mmatrix + 1.0_wp/(2.0_wp*timestep)*cmatrix    	
         	case(.true.)
                transmatrix = mmatrix*1.0_wp/(timestep**2.0_wp)     
        end select
        

       	!Enforce BCs and factorize left side(transmatrix)
        call enforce(transmatrix)
        call factor(transmatrix)   

		MD = 0
        CD = 0
		! MAIN EXPLICIT DIRECT INTEGRATION LOOP
        do n = 1, no_timestep-1
!$$$$$$             print *, ""
!$$$$$$             print *, "Time step number: ",n
!$$$$$$             print *, ""
            ! Save old variables
            MD_old = MD
            CD_old = CD
			
			! Build P (if transient force is used)
            if ((transient_force) .and. (n == 2)) then
                p = 0
            end if
    
            ! Build MD (multiplied elementwise)
            select case(lumped_mass)
            	case(.true.)
                	call mmul(D_all(:,n),MD,3) !MD
                case(.false.)
					call mmul(D_all(:,n),MD,2) !MD
            end select
            ! Build KD(multiplied elementwise)
            call mmul(D_all(:,n),KD,1)
            ! Build CD(multiplied elementwise)
            select case (boundary_c)
                case(.false.)
                    CD = alpha*MD + beta*KD
                case(.true.)
                	call mmul(D_all(:,n),CD,4)
        	end select

            ! Build right side (put directly into D_all(:,n+1) for solving)
            select case(half_Step)
              	case(.false.)
					D_all(:,n+1) = P - KD + 2.0_wp/(timestep**2)*MD - MD_old/(timestep**2) + CD_old/(2*timestep) ! (11.12-3)
                case(.true.)
                	D_all(:,n+1) = P - KD + 2.0_wp/(timestep**2)*MD - CD/timestep - MD_old/(timestep**2) + CD_old/(timestep) !(11.12-8)
			end select
			
			! Enforce BC for højresiden
            do i = 1, nb
                idof = int(2*(bound(i,1)-1) + bound(i,2))
            	D_all(idof,n+1) = bound(i, 3)
            end do
            
			! Løs ligningssystemet
            call solve(transmatrix,D_all(:,n+1))

			! Plot
            d(1:neqn) = D_all(1:neqn,n+1)
            call plot(deformed, wait=.false.)
        																							
		end do

!$$$$$$         do n = 1, no_timestep
!$$$$$$             print *, "---------------------------------------------------------------------"
!$$$$$$             print *, "Displacement no:", n
!$$$$$$             print *,D_all(:,n)
!$$$$$$         end do

      	!Print Y-displacement
!$$$$$$         print *, d(146*2-1)
!$$$$$$         do i = 1,no_timestep-1
!$$$$$$           print *, D_all(146*2-1,i)
!$$$$$$         end do

		!Output results
        if (matlab_output) then
            write(t,'(f20.4)') timestep
            print '(A)', trim(t)
            outfile_path = trim(filename)//'_t'//trim(t)//'.m'
            print *, trim(outfile_path)
            open(10, file=trim(outfile_path))
            write(10,*) 'clear all;close all;F=false;T=true'
    
            write(10,*) 'alpha =',alpha
            write(10,*) 'beta =',beta
            write(10,*) 'timestep =',timestep
            write(10,*) 'no_timestep =',no_timestep
            write(10,*) 'lumped_mass =',lumped_mass
            write(10,*) 'half_step =',half_step
            write(10,*) 'displ = ['
            do i = 1,no_timestep
                write(10,*) D_all(364*2-1,i)
            end do
            write(10,*) '];'
            write(10,*) 'plot(0:timestep:timestep*(no_timestep-1),displ)'
            close(10)
 		end if
     end subroutine trans
 !
 !--------------------------------------------------------------------------------------------------
 !
     subroutine eigen
 
         !! This subroutine calculates the eigenvalue and eigenfrequencies
 
         use fedata
         use numeth
         use processor
         use plane42
         
         integer :: e, pp, pmax, n, nen, i, idof, q, j, m, mtype
         !real(wp), dimension(:), allocatable :: plotval
         real(wp):: Xnew(neqn), Ynew(neqn), Xold(neqn), Yold(neqn), r, z(neqn,no_eigenvalues),  lambda(no_eigenvalues)
 
         real(wp), dimension(mdim, mdim) :: me 
         real(wp), dimension(mdim) :: xe
         integer, dimension(mdim) :: edof
         real(wp) :: young!, area
         real(wp) :: nu, dens, thk, c
         
         lambda = 0.0_wp
         mtype = 3
 
         
         if (.not. isoparametric) then
             print *, "nodal analysis is not implemented for non-isoparametric elements"
             stop
         end if
 
 
         ! Build load-vector
         call buildload
 
         ! Build stiffness matrix
         call buildstiff
         
         ! Remove rigid body modes
         call enforce(kmatrix)
 
         if (.not. banded) then
             call factor(kmatrix)
         else    
             call bfactor(kmatrix)
         end if
 
      do m = 1, no_eigenvalues
         ! initial eigenvector
         Xnew = 1
         
         ! Calculate Y0 from X0 and M (Element multiplication)
         call mmul(Xnew,Ynew,mtype)
         
 
         !compute Z to find higer natural frequencies, only after first run will this loop activate
         do j=1, m-1
           call mmul(D_eig(:,j),z(:,j),mtype)
         end do
 
 
         ! iterate start
         pmax = 5000
         do pp = 1, pmax
           Xold = Xnew
           Yold = Ynew
           
           !enforce boundary conditions on Y as if on p in displacement calculations
           if (.not. banded) then
               do i = 1, nb !iterate over all bcs
                   idof = int(2*(bound(i,1)-1) + bound(i,2)) !find dof of bc
                   Ynew(idof) = 0.0_wp
               end do
           else
               if (.not. penalty) then
                   do i = 1, nb !iterate over all bcs
                       idof = int(2*(bound(i,1)-1) + bound(i,2)) !find dof of bc
                       Ynew(idof) = 0.0_wp
                   end do
               else
                    print *, "Penalty method not implemented for banded stiffness matrix fea.f90(enforce subroutine)"
               end if
           end if
           
           !solve for X according to banded or not banded
           Xnew(1:neqn) = Yold(1:neqn)
           if (.not. banded) then
             call solve(kmatrix, Xnew)
           else    
             call bsolve(kmatrix, Xnew)
           end if
           
           !Enforce boundary condition on Xnew
           do i = 1, nb !iterate over all bcs
             idof = int(2*(bound(i,1)-1) + bound(i,2)) !find dof of bc
             Xnew(idof) = 0.0_wp
           end do
           
           !compute c used for orthogonalize
           do j=1, m-1
             c = dot_product(Xnew,Z(:,j))
 !$$$$$$             print *, "c: ", c
             Xnew = Xnew - c*D_eig(:,j)
           end do
             
           !Calculate Y from X(eigenvector) found in solve
           call mmul(Xnew,Ynew,mtype)
           
           ! Calculate r
           r = sqrt(dot_product(Xnew,Ynew))
           
           ! Calculate Y for this iteration
           Ynew = Ynew / r
           
           ! Check if precision for inverse iteration is met
           !print *, "precision: ", sqrt(dot_product(Xnew-Xold,Xnew-Xold))/sqrt(dot_product(Xnew,Xnew))
           if (sqrt(dot_product(Xnew-Xold,Xnew-Xold))/sqrt(dot_product(Xnew,Xnew)) < estop) then
             print *, "Number of iterations:", pp
             
             exit
           end if
           
           ! Inform if precision was not met
           if (pp == pmax) then
             print *, "Max iterations have been met. fea.f90(eigen subroutine)"
           end if     
         end do
         
         ! eigenvector
         D_eig(:,m) = Xnew/r
 
         ! eigenfreq
         lambda(m) = dot_product(Xnew,Yold)/(r**2)
         print *, "Natural frequencies number:",m , sqrt(lambda(m))
         
         ! Plot deformed shape
         !call plot(eigenmode, eigenfreq=lambda(m), eigenvector=D_eig(:,m))
         
         !! plot data to matlabfile
         !call plot(eigenmode, eigenfreq=lambda(m), eigenvector=0.5*D_eig(:,m), tend=0.1_wp,tdelta=10.0_wp, &
         !device=Matlab, title='Eigenmode')
         
         !!Plot with time control
         !call plot(eigenmode, eigenfreq=lambda(m), eigenvector=D_eig(:,m), tend=10.0_wp,tdelta=0.01_wp)
      end do
     end subroutine eigen
 
 !
 !--------------------------------------------------------------------------------------------------
 !
 
     subroutine mmul(invector,outvector, mtype)
 
     use fedata
     use plane42
     
     integer, intent(in) :: mtype
     real(wp), intent(in) :: invector(neqn)
     real(wp), intent(out) :: outvector(neqn)
     
     integer :: e, n, nen
     integer, dimension(mdim) :: edof
     real(wp) :: A, s, dens, thk, young, nu, xe(mdim), e_mat(mdim,mdim), e_vect(mdim)
 
     outvector(:) = 0.0_wp
     
     do e = 1, ne
             ! Find coordinates
             nen = element(e)%numnode
             do n = 1, nen
                  xe(2*n-1) = x(element(e)%ix(n),1)
                  xe(2*n  ) = x(element(e)%ix(n),2)
                  edof(2*n-1) = 2 * element(e)%ix(n) - 1  
                  edof(2*n)   = 2 * element(e)%ix(n)
             end do
 
             ! Load element properties
             thk = mprop(element(e)%mat)%thk ! thk - element thickness
             dens = mprop(element(e)%mat)%dens ! dens - element Density
             young = mprop(element(e)%mat)%young ! young - E/young's modulus
             nu = mprop(element(e)%mat)%nu ! nu - Poisson's ratio
             
             ! Multiply matrix according to chosen mtype
             select case (mtype)
                 case(1)
                     ! Multiply stiffness matrix with "invect" elementwise
                     call plane42_ke(xe, young, nu, thk, e_mat)
                     outvector(edof) = outvector(edof) + matmul(e_mat,invector(edof))
                 case(2)
                     ! Multiply consistent mass matrix with "invect" elementwise
                     call plane42_me(xe, thk, dens, e_mat)
                     outvector(edof) = outvector(edof) + matmul(e_mat,invector(edof))
                 case(3)
                     ! Multiply lumped mass matrix with "invect" elementwise
                     call plane42_lme(xe, thk, dens, e_vect)
                     do n = 1,size(edof)
                         outvector(edof(n)) = outvector(edof(n)) + e_vect(n)*invector(edof(n))
                     end do
 				case(4)
                	! Multiply damping for BCs with "invect" elementwise
                    call plane42_cebc(e, e_mat)
                    outvector(edof) = outvector(edof) + matmul(e_mat,invector(edof))
             end select
         end do
 
     end subroutine mmul
 
 !
 !--------------------------------------------------------------------------------------------------
 !
     subroutine buildload
 
         !! This subroutine builds the global load vector
 
         use fedata
         use plane42
         use plane42rect
 
         integer :: i, e, nen, pdof, n, eface
         ! Hint for continuum elements:
         integer, dimension(mdim) :: edof
         real(wp), dimension(mdim) :: xe
         real(wp), dimension(mdim) :: re
         real(wp) :: thk, fe
 
         ! Build load vector
         p(1:neqn) = 0
         do i = 1, np
             ! print load for debugging
             ! print *, 'Processing load:', i, '   Case: ', int(loads(i, 1))
             
             select case(int(loads(i, 1)))
             case( 1 )
                 ! Build nodal load contribution
                 pdof = loads(i, 2)*2 - 2 + loads(i,3)
                 p(pdof) = p(pdof) + loads(i,4)            
             case( 2 )
             	 if (.not. antype == 'trans') then
                     ! Build uniformly distributed surface (pressure) load contribution
                     e = loads(i,2)
                     ! Find coordinates
                     nen = element(e)%numnode
                     do n = 1, nen
                          xe(2*n-1) = x(element(e)%ix(n),1)
                          xe(2*n  ) = x(element(e)%ix(n),2)
                          edof(2*n-1) = 2 * element(e)%ix(n) - 1  
                          edof(2*n)   = 2 * element(e)%ix(n)
                          ! print *, 'node',n,'x',edof(2*n-1)
                          ! print *, 'node',n,'y',edof(2*n)
                     end do
                     eface = loads(i,3)
                     fe = loads(i,4)
                     thk = mprop(element(e)%mat)%thk ! thk - element thickness
                     select case (isoparametric)
                       case(.true.)
                         call plane42_re(xe, eface, fe, thk, re) ! returns re - Element load vector (output)
                       case(.false.)
                         call plane42rect_re(xe, eface, fe, thk, re) ! returns re - Element load vector (output)
                     end select
                     p(edof) = p(edof) + re
                     print *, "ithink antype is not trans"
                 end if
                 
                 
             case default
                 print *, 'ERROR in fea/buildload'
                 print *, 'Load type not known'
                 stop
             end select
         end do
     end subroutine buildload
 !
 !--------------------------------------------------------------------------------------------------
 !
     subroutine buildstiff
 
         !! This subroutine builds the global stiffness matrix from
         !! the local element stiffness matrices
 
         use fedata
         use link1
         use plane42
         use plane42rect
 
         integer :: e, i, j, ibw, jbw
         integer :: nen
 ! Hint for system matrix in band form:
 !        integer :: irow, icol
         integer, dimension(mdim) :: edof
         real(wp), dimension(mdim) :: xe
         real(wp), dimension(mdim, mdim) :: ke 
         real(wp) :: young, area, thk, nu
 
         ! Reset stiffness matrix
         kmatrix = 0
 
         do e = 1, ne
           
             ! Find coordinates and degrees of freedom
             !nen = element(e)%numnode
             do i = 1, element(e)%numnode
                  xe(2*i-1) = x(element(e)%ix(i),1)
                  xe(2*i  ) = x(element(e)%ix(i),2)
                  edof(2*i-1) = 2 * element(e)%ix(i) - 1  
                  edof(2*i)   = 2 * element(e)%ix(i)
             end do
 
             ! Gather material properties and find element stiffness matrix
             select case( element(e)%id )
             case( 1 ) ! element is of type truss
                  young = mprop(element(e)%mat)%young
                  area  = mprop(element(e)%mat)%area
                  call link1_ke(xe, young, area, ke)
             case( 2 ) ! element is of type 4-node
                  ! xe - Nodal coordinates of this element in undeformed configuration
                  young = mprop(element(e)%mat)%young ! young - E/young's modulus
                  nu = mprop(element(e)%mat)%nu ! nu - Poisson's ratio
                  thk = mprop(element(e)%mat)%thk ! thk - element thickness
                  select case (isoparametric)
                      case(.true.)
                        call plane42_ke(xe, young, nu, thk, ke) ! returns ke - Stiffness matrix
                      case(.false.)
                        call plane42rect_ke(xe, young, nu, thk, ke) ! returns re - Element load vector (output)
                 end select
                  
             end select
 
             ! Assemble into global matrix
             if (.not. banded) then
                   kmatrix(edof,edof) = kmatrix(edof,edof) + ke
             else
 !PRINT
 !$$$$$$                   print *, '-------- Element:',e,'--------'
 !$$$$$$                   print *, 'edofs: ',edof
                   do i = 1, 2*element(e)%numnode
                       do j = 1, 2*element(e)%numnode
                           ! ibw and jbw converts the indices of dofs to indices in the banded matrix
                           if ((edof(i)-edof(j)) >= 0) then
                               ibw = edof(i) - edof(j) + 1
                               jbw = min(edof(j),edof(i))
                               !print *, "ij:",i,j," edof(i)edof(j)",edof(i),edof(j)," ibwjbw:",ibw,jbw
                               kmatrix(ibw,jbw) = kmatrix(ibw,jbw) + ke(i,j)
                           end if
                       end do
                   end do
             end if 
         end do
 
 !$$$$$$ !        print k
 !$$$$$$          do i = 1,bw
 !$$$$$$            print "(24(f4.2,tr1))",kmatrix(i,1:neqn)
 !$$$$$$          end do
 !$$$$$$         pause
     end subroutine buildstiff
 !
 !--------------------------------------------------------------------------------------------------
 !
     subroutine enforce(enforce_matrix)
 
         !! This subroutine enforces the support boundary conditions
 
         use fedata
 
         integer :: i, j, idof
         real(wp) :: penal
         real(wp), intent(out) :: enforce_matrix(neqn,neqn)
 
         ! Correct for supports
         if (.not. banded) then
             if (.not. penalty) then
                 do i = 1, nb
                     idof = int(2*(bound(i,1)-1) + bound(i,2))
                     p(1:neqn) = p(1:neqn) - enforce_matrix(1:neqn, idof) * bound(i, 3)
                     p(idof) = bound(i, 3)
                     enforce_matrix(1:neqn, idof) = 0
                     enforce_matrix(idof, 1:neqn) = 0
                     enforce_matrix(idof, idof) = 1
                 end do
             else
                 penal = penalty_fac*maxval(enforce_matrix)
                 do i = 1, nb
                     idof = int(2*(bound(i,1)-1) + bound(i,2))
                     enforce_matrix(idof, idof) = enforce_matrix(idof, idof) + penal
                     p(idof) = penal * bound(i, 3)  
                 end do 
             end if
         else
           if (.not. penalty) then
             do i = 1, nb
                 ! loop igennem boundary conditions
                 idof = int(2*(bound(i,1)-1) + bound(i,2))
                 !p(1:neqn) = p(1:neqn) - enforce_matrix(1:neqn, idof) * bound(i, 3)
                 p(idof) = bound(i, 3)
 
                 ! Sæt 1 i (idof,idof)
                 enforce_matrix(1, idof) = 1
 
                 ! Sæt 0 i banded k svarende til kolonnen i gammel enforce_matrix
                 enforce_matrix(2:bw, idof) = 0
                 
                 ! Sæt 0 i banded k svarende til rækken i gammel enforce_matrix
                 do j = 2,bw
 !PRINT              print *, 'index',j,idof-(j-1)
                     if (idof-(j-1) > 0) then
                         enforce_matrix(j, idof-(j-1)) = 0
                     else
                         exit
                     end if
                 end do
              end do
            else
              print *, "Penalty method not implemented for banded stiffness matrix fea.f90(enforce subroutine)"
            end if
         end if
     end subroutine enforce
 !
 !--------------------------------------------------------------------------------------------------
 !
     subroutine recover
 
         !! This subroutine recovers the element stress, element strain, 
         !! and nodal reaction forces
 
         use fedata
         use link1
         use plane42
         use plane42rect
 
         integer :: e, i, nen
         integer :: edof(mdim)
         real(wp), dimension(mdim) :: xe, de
         real(wp), dimension(mdim, mdim) :: ke
         real(wp) :: young, area 
         real(wp):: nu, dens, thk
         real(wp), dimension(3) :: estrain, estress
         real(wp):: stress_d, stress_1, stress_2, evm_stress, s, c
 
         ! PRINT Global Compliance
         print *, "Global compliance", dot_product(d,p)
 
         ! Reset force vector
         p = 0
 
         do e = 1, ne
           
             ! Find coordinates etc...
             nen = element(e)%numnode
             do i = 1,nen
                 xe(2*i-1) = x(element(e)%ix(i), 1)
                 xe(2*i)   = x(element(e)%ix(i), 2)
                 edof(2*i-1) = 2 * element(e)%ix(i) - 1
                 edof(2*i)   = 2 * element(e)%ix(i)
                 de(2*i-1) = d(edof(2*i-1))
                 de(2*i)   = d(edof(2*i))
             end do
 
             ! Find stress and strain
             select case( element(e)%id )
             case( 1 )
                 young = mprop(element(e)%mat)%young
                 area  = mprop(element(e)%mat)%area
                 call link1_ke(xe, young, area, ke)
                 p(edof(1:2*nen)) = p(edof(1:2*nen)) + matmul(ke(1:2*nen,1:2*nen), de(1:2*nen))
                 call link1_ss(xe, de, young, estress, estrain)
                 stress(e, 1:3) = estress
                 strain(e, 1:3) = estrain
             case( 2 )
                 young = mprop(element(e)%mat)%young
                 nu = mprop(element(e)%mat)%nu
                 select case (isoparametric)
                   case(.true.)
                     call plane42_ss(xe, de, young, nu, estress, estrain)
                   case(.false.)
                     call plane42rect_ss(xe, de, young, nu, estress, estrain, element_evaluation_x, element_evaluation_y)
                 end select
                 
                 stress(e, 1:3) = estress
                 strain(e, 1:3) = estrain
                 
                 stress_d = sqrt((estress(1)-estress(2)/0.5)**2+estress(3)**2)
                 stress_1 = 0.5*(estress(1)+estress(2))+ stress_d
                 stress_2 = 0.5*(estress(1)+estress(2))- stress_d   
                 evm_stress = sqrt(stress_1**2+stress_2**2-stress_1*stress_2)
                 
             
                 vm_stress(e) = evm_stress   
 
 
 
                 s = -2*estress(3) / (stress_1-stress_2)
                 c = (estress(1)-estress(2)) / (stress_1-stress_2)
                 psi(e) = atan2(s,c)/2
                 
                 
 !$$$$$$                 if (sqrt(x(element(e)%ix(1),1)**2+x(element(e)%ix(1),2)**2) < 1.201) then
 !$$$$$$                   !print *, "distance:",sqrt(x(element(e)%ix(2),1)**2+x(element(e)%ix(2),2)**2)
 !$$$$$$                   print *, "X:", x(element(e)%ix(2),1),"Y:",x(element(e)%ix(2),2)
 !$$$$$$                   !print *, 'stress vm: ',stress_vm,"direction:",psi
 !$$$$$$                 end if
             end select
             
 ! PRINT Point A and B stresses
             if(x(element(e)%ix(3),1)==3.5 .AND. x(element(e)%ix(3),2)==0.0) then
                 print *, "Point A stress", stress(e, 1:3)
                 print *, "Point A Von mises stress", evm_stress
             end if
 !$$$$$$             if(x(element(e)%ix(4),1)==1.0 .AND. x(element(e)%ix(4),2)==1.0) then
 !$$$$$$                 print *, "Point B stress", stress(e, 1:3)
 !$$$$$$             end if
         end do
         
     end subroutine recover
 end module fea
 