module plane42 

    !! This module contains subroutines specific to the plane42 element 
    !!        
    !! The plane42 element has 4 nodes. Each node has 2 degrees-of-freedom,
    !! namely, displacement along the \(x\)- and \(y\)-coordinate directions.
    !!
    !! The nodes are numbered counter-clockwise as indicated below. The
    !! figure also shows how the element edges are labelled. For example,
    !! edge 1 is the edge between element node 1 and 2.
    !!
    !!       N4    E3    N3
    !!          o------o 
    !!          |      |
    !!       E4 |      | E2
    !!          |      |
    !!          o------o
    !!       N1    E1    N2
    !!             
    !!
    !! `N1` = element node 1, `N2` = element node 2, etc  
    !! `E1` = element face 1, `E2` = element face 2, etc  
    
    use types
    implicit none
    save
    
    private
    public :: plane42_N, plane42_shape, plane42_me, plane42_ke, plane42_re, plane42_ss, plane42_lme, plane42_cebc
contains

	subroutine plane42_N(xi, eta, N)
        !! This subroutine calculates N for one gauss point in a plane42 element

        use fedata
        
        real(wp), intent(in) :: xi, eta
        real(wp), intent(out) :: N(2,8)
        real(wp) :: L(3,4)
        integer :: i
            
        ! Calculate L
        L = 0
        L(1,1) = 1
        L(2,4) = 1
        L(3,2) = 1
        L(3,3) = 1

        ! Calculate N
        N = 0.0_wp
        N(1,1) = 0.25_wp*(1-eta)*(1-xi)
        N(2,2) = 0.25_wp*(1-eta)*(1-xi)
        N(1,3) = 0.25_wp*(1+eta)*(1-xi)
        N(2,4) = 0.25_wp*(1+eta)*(1-xi)
        N(1,5) = 0.25_wp*(1+eta)*(1+xi)
        N(2,6) = 0.25_wp*(1+eta)*(1+xi)
        N(1,7) = 0.25_wp*(1-eta)*(1+xi)
        N(2,8) = 0.25_wp*(1-eta)*(1+xi)

    end subroutine plane42_N

!----------------------------------------------------------------------------------------------
        
	subroutine plane42_shape(xi, eta, N, B, det_J, xe)
        !! This subroutine calculates N, B, Gamma and |J| for one gauss point 
        !! in a plane42 element

        use fedata
        
		real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42]]
        real(wp), intent(in) :: xi, eta
        	!! coordinates for gauss point
        real(wp), intent(out) :: B(3,8), det_J, N(2,8)
			!!
        real(wp) :: L(3,4), Gamma(2,2), Gammatilde(4,4), Ntilde(4,8)
        	!!
        real(wp) :: N1dxi,N2dxi,N3dxi,N4dxi, &
        			N1deta,N2deta,N3deta,N4deta
        real(wp) :: J11, J12, J21, J22
        
        integer :: i
            
        ! Calculate L
        L = 0
        L(1,1) = 1
        L(2,4) = 1
        L(3,2) = 1
        L(3,3) = 1

        ! Calculate N
        N = 0.0_wp
        N(1,1) = 0.25_wp*(1-eta)*(1-xi)
        N(2,2) = 0.25_wp*(1-eta)*(1-xi)
        N(1,3) = 0.25_wp*(1+eta)*(1-xi)
        N(2,4) = 0.25_wp*(1+eta)*(1-xi)
        N(1,5) = 0.25_wp*(1+eta)*(1+xi)
        N(2,6) = 0.25_wp*(1+eta)*(1+xi)
        N(1,7) = 0.25_wp*(1-eta)*(1+xi)
        N(2,8) = 0.25_wp*(1-eta)*(1+xi)
        
        ! Calculate Ntilde
		N1dxi = -0.25_wp+0.25_wp*eta
        N1deta = -0.25_wp+0.25_wp*xi
        N2dxi =  0.25_wp-0.25_wp*eta
        N2deta = -0.25_wp-0.25_wp*xi
        N3dxi =  0.25_wp+0.25_wp*eta
        N3deta =  0.25_wp+0.25_wp*xi
        N4dxi = -0.25_wp-0.25_wp*eta
        N4deta =  0.25_wp-0.25_wp*xi

        Ntilde = 0
        Ntilde(1,1) = N1dxi
        Ntilde(2,1) = N1deta
        Ntilde(3,2) = N1dxi
        Ntilde(4,2) = N1deta
       
        Ntilde(1,3) = N2dxi
        Ntilde(2,3) = N2deta
        Ntilde(3,4) = N2dxi
        Ntilde(4,4) = N2deta

        Ntilde(1,5) = N3dxi
        Ntilde(2,5) = N3deta
        Ntilde(3,6) = N3dxi
        Ntilde(4,6) = N3deta

        Ntilde(1,7) = N4dxi
        Ntilde(2,7) = N4deta
        Ntilde(3,8) = N4dxi
        Ntilde(4,8) = N4deta

        ! Calculate J and |J|
        J11 = dot_product(Ntilde(1,:),xe)
        J21 = dot_product(Ntilde(2,:),xe)
        J12 = dot_product(Ntilde(3,:),xe)
        J22 = dot_product(Ntilde(4,:),xe)
        det_J = J11*J22-J12*J21
!$$$$$$         print *, '----------------------------------'
!$$$$$$         print *, 'J:'
!$$$$$$         print '(24(f4.2,tr1))', J11,J12
!$$$$$$         print '(24(f4.2,tr1))', J21,J22
!$$$$$$         print *, 'det_J',det_J

        ! Calculate Gammatilde
        Gamma = 0
        Gamma(1,1) = 1/det_J*J22
        Gamma(1,2) = 1/det_J*(-J12)
        Gamma(2,1) = 1/det_J*(-J21)
        Gamma(2,2) = 1/det_J*J11

        Gammatilde = 0
        Gammatilde(1,1) = Gamma(1,1)
        Gammatilde(1,2) = Gamma(1,2)
        Gammatilde(2,1) = Gamma(2,1)
        Gammatilde(2,2) = Gamma(2,2)
        Gammatilde(3,3) = Gamma(1,1)
        Gammatilde(3,4) = Gamma(1,2)
        Gammatilde(4,3) = Gamma(2,1)
        Gammatilde(4,4) = Gamma(2,2)
!PRINT      
!$$$$$$         print *, 'Gammatilde:'
!$$$$$$         do i = 1,4
!$$$$$$             print "(24(f4.2,tr1))",Gammatilde(i,:)
!$$$$$$         end do
        
        ! Calculate B
		B = matmul(matmul(L,Gammatilde),Ntilde)
!PRINT
!$$$$$$         print *, 'B',shape(B)
!$$$$$$         do i = 1,4
!$$$$$$             print "(24(f4.2,tr1))",B(i,:)
!$$$$$$         end do
        

    end subroutine plane42_shape

!----------------------------------------------------------------------------------------------

    subroutine plane42_lme(xe, thk, dens, me)

        !! This subroutine constructs the mass matrix for
        !! a rectangular 4-noded quad element.

        use fedata

        
        real(wp), intent(in) :: thk, dens
        real(wp), dimension(:), intent(in) :: xe
        real(wp), dimension(:), intent(out) :: me
        integer :: i, j, ike
		real(wp) :: loc(3), w(3), xi, eta
        real(wp) :: NtN(8), N(2,8), det_J, A
        real(wp) :: Ntilde(4,8)
        real(wp) :: N1dxi,N2dxi,N3dxi,N4dxi, &
        			N1deta,N2deta,N3deta,N4deta
        real(wp) :: J11, J12, J21, J22
        
		! Reset variables
		loc = 0
        w = 0
        me = 0
        
		! Define variable Gauss points and weighting
        select case (no_gauss_points_k)
            case (1)
            	loc(1) = 0
                w(1) = 2
            case(2)
                loc(1) = 1_wp/sqrt(3.0_wp)
                loc(2) = -1_wp/sqrt(3.0_wp)
                w(1) = 1
                w(2) = 1
            case(3)
            	loc(1) = sqrt(0.6_wp)
                loc(2) = 0
                loc(3) = -sqrt(0.6_wp)
                w(1) = 5.0_wp/9.0_wp
                w(2) = 8.0_wp/9.0_wp
                w(3) = -w(1)
            end select

        ! Calculate me for each gauss point
        do i = 1, no_gauss_points_k
            do j = 1, no_gauss_points_k
                xi = loc(i)
                eta = loc(j)
        
                ! Calculate N^T*N
                NtN = 0.0_wp
                NtN(1) = (0.25_wp*(1-eta)*(1-xi))**2
                NtN(2) = (0.25_wp*(1-eta)*(1-xi))**2
                NtN(3) = (0.25_wp*(1+eta)*(1-xi))**2
                NtN(4) = (0.25_wp*(1+eta)*(1-xi))**2
                NtN(5) = (0.25_wp*(1+eta)*(1+xi))**2
                NtN(6) = (0.25_wp*(1+eta)*(1+xi))**2
                NtN(7) = (0.25_wp*(1-eta)*(1+xi))**2
                NtN(8) = (0.25_wp*(1-eta)*(1+xi))**2
        		! Sum Numerical integration of me

                ! Calculate Ntilde
                N1dxi = -0.25_wp+0.25_wp*eta
                N1deta = -0.25_wp+0.25_wp*xi
                N2dxi =  0.25_wp-0.25_wp*eta
                N2deta = -0.25_wp-0.25_wp*xi
                N3dxi =  0.25_wp+0.25_wp*eta
                N3deta =  0.25_wp+0.25_wp*xi
                N4dxi = -0.25_wp-0.25_wp*eta
                N4deta =  0.25_wp-0.25_wp*xi
        
                Ntilde = 0
                Ntilde(1,1) = N1dxi
                Ntilde(2,1) = N1deta
                Ntilde(3,2) = N1dxi
                Ntilde(4,2) = N1deta
               
                Ntilde(1,3) = N2dxi
                Ntilde(2,3) = N2deta
                Ntilde(3,4) = N2dxi
                Ntilde(4,4) = N2deta
        
                Ntilde(1,5) = N3dxi
                Ntilde(2,5) = N3deta
                Ntilde(3,6) = N3dxi
                Ntilde(4,6) = N3deta
        
                Ntilde(1,7) = N4dxi
                Ntilde(2,7) = N4deta
                Ntilde(3,8) = N4dxi
                Ntilde(4,8) = N4deta
        
                ! Calculate J and |J|
                J11 = dot_product(Ntilde(1,:),xe)
                J21 = dot_product(Ntilde(2,:),xe)
                J12 = dot_product(Ntilde(3,:),xe)
                J22 = dot_product(Ntilde(4,:),xe)
                det_J = J11*J22-J12*J21
                me = me + w(i)*w(j)*thk*dens*NtN*det_J
        	end do
        end do
        ! Scale lumped mass matrix to correct mass
        A = det_J*4
        
        me = me * dens*thk*A/sum(me)
        
    end subroutine plane42_lme

!----------------------------------------------------------------------------------------------

    subroutine plane42_me(xe, thk, dens, me)

        !! This subroutine constructs the mass matrix for
        !! a rectangular 4-noded quad element.

        use fedata

        
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), intent(in) :: dens
            !! Density of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42]]
        real(wp), dimension(:,:), intent(out) :: me
            !! Stiffness matrix
        integer :: i, j, ike
		real(wp) :: loc(3), w(3), xi, eta
        	!! Gauss numerical integration
        real(wp), dimension(3, mdim) :: B
        real(wp) :: det_J
        real(wp) :: N(2,8)
        	!! Outputs of shape subroutime
        
		! Reset variables
		loc = 0
        w = 0
        me = 0
        
		! Define variable Gauss points and weighting
        select case (no_gauss_points_k)
            case (1)
            	loc(1) = 0
                w(1) = 2
            case(2)
                loc(1) = 1_wp/sqrt(3.0_wp)
                loc(2) = -1_wp/sqrt(3.0_wp)
                w(1) = 1
                w(2) = 1
            case(3)
            	loc(1) = sqrt(0.6_wp)
                loc(2) = 0
                loc(3) = -sqrt(0.6_wp)
                w(1) = 5.0_wp/9.0_wp
                w(2) = 8.0_wp/9.0_wp
                w(3) = -w(1)
            end select

        ! Calculate me for each gauss point
        do i = 1, no_gauss_points_k
            do j = 1, no_gauss_points_k
                xi = loc(i)
                eta = loc(j)
                ! Calculate N, B and detJ with shape subroutine
				call plane42_shape(xi, eta, N, B, det_J, xe)
        		! Sum Numerical integration of me
                me = me + w(i)*w(j)*thk*dens*matmul(transpose(N),N)*det_J
        	end do
        end do
    end subroutine plane42_me

!----------------------------------------------------------------------------------------------

    subroutine plane42_ke(xe, young, nu, thk, ke)

        !! This subroutine constructs the stiffness matrix for
        !! a rectangular 4-noded quad element.

        use fedata

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu
            !! Poisson's Ratio for this element
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration
            !!
            !! * `xe(1:2)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(3:4)` = \((x,y)\)-coordinates of element node 1
            !! * `xe(5:6)` = \((x,y)\)-coordinates of element node 2
            !! * `xe(7:8)` = \((x,y)\)-coordinates of element node 2
            !!
            !! See also [[plane42]]
        real(wp), dimension(:,:), intent(out) :: ke
            !! Stiffness matrix
        real(wp) :: cmat(3,3), fact
			!! C, constitutive matrix
        integer :: i, j, ike
		real(wp) :: loc(3), w(3), xi, eta
        	!! Gauss numerical integration
        real(wp), dimension(3, mdim) :: B
        real(wp) :: det_J
        real(wp) :: N(2,8)
        	!! Outputs of shape subroutime
        
		! Reset variables
		loc = 0
        w = 0
        ke = 0
        
		! Define variable Gauss points and weighting
        select case (no_gauss_points_k)
            case (1)
            	loc(1) = 0
                w(1) = 2
            case(2)
                loc(1) = 1_wp/sqrt(3.0_wp)
                loc(2) = -1_wp/sqrt(3.0_wp)
                w(1) = 1
                w(2) = 1
            case(3)
            	loc(1) = sqrt(0.6_wp)
                loc(2) = 0
                loc(3) = -sqrt(0.6_wp)
                w(1) = 5.0_wp/9.0_wp
                w(2) = 8.0_wp/9.0_wp
                w(3) = -w(1)
            end select

        ! build constitutive matrix (plane stress)
        cmat = 0
        fact = young/(1-nu**2)
        cmat(1,1) = fact
        cmat(1,2) = fact*nu
        cmat(2,1) = fact*nu
        cmat(2,2) = fact
        cmat(3,3) = fact*(1-nu)/2

        ! Calculate ke for each gauss point
        do i = 1, no_gauss_points_k
            do j = 1, no_gauss_points_k
                xi = loc(i)
                eta = loc(j)
                ! Calculate N, B and detJ with shape subroutine
				call plane42_shape(xi, eta, N, B, det_J, xe)
        		! Calculate K
                ke = ke + w(i)*w(j)*thk*matmul(matmul(transpose(B),cmat),B)*det_J
        	end do
        end do
!PRINT
!$$$$$$         print *, 'ke:'
!$$$$$$         do ike = 1,8
!$$$$$$             print "(24(f4.2,tr1))",ke(ike,:)
!$$$$$$         end do
        
    end subroutine plane42_ke
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_re(xe, eface, fe, thk, re)

        !! This subroutine computes the element load vector due
        !! to surface traction (traction is always perpendicular
        !! to element face).

        ! inputs:
        ! xe		nodal coordniates
        ! eface		face number
        ! fe		pressure value
        ! thk		thickness
        ! re		element load vector (output)
        

        integer, intent(in) :: eface
            !! Element face where traction (pressure) is applied

        real(wp), intent(in) :: fe
            !! Value of surface traction (pressure)
        real(wp), intent(in) :: thk
            !! Thickness of this element
        real(wp), dimension(:), intent(in) :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42_ke]])
        real(wp), intent(out) :: re(8)
            !! Element force vector
            !!
            !! * `re(1:2)` = \((f_x^1, f_y^1)\) force at element node 1 in \(x\)- and y-direction
            !! * `re(3:4)` = \((f_x^2, f_y^2)\) force at element node 1 in \(x\)- and y-direction
            !! * etc...
        real(wp) :: aa, bb, nface(2,8), f(2)

        aa = (xe(3)-xe(1))/2
        bb = (xe(8)-xe(2))/2

        nface = 0
        f = 0
        if (eface == 1) then
            nface(1,1) = aa
            nface(1,3) = aa
            nface(2,2) = aa
            nface(2,4) = aa
            f(2) = fe
        elseif (eface == 2) then
            nface(1,3) = bb
            nface(1,5) = bb
            nface(2,4) = bb
            nface(2,6) = bb
            f(1) = -fe
        elseif (eface == 3) then
            nface(1,5) = aa
            nface(1,7) = aa
            nface(2,6) = aa
            nface(2,8) = aa
            f(2) = -fe
        elseif (eface == 4) then
            nface(1,1) = bb
            nface(1,7) = bb
            nface(2,2) = bb
            nface(2,8) = bb
            f(1) = fe
        endif
        
        re = matmul(transpose(nface), f) * thk
        
    end subroutine plane42_re
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_cebc(e, cebc)

        !! This subroutine computes the element damping vector due
        !! to absorbing BCs (modeled in flextract as surface tractions)
        !! (traction is always perpendicular to element face).

        use fedata
        use fortranhacks

        integer, intent(in) :: e
        real(wp), intent(out) :: cebc(8,8)
        integer :: plane_stress_or_strain, eface, n, i, nen
        real(wp) :: thk, young, nu, dens
        real(wp) :: xe(mdim), z(2,2), aa, bb, nface(2,8), cp, cs, surface_normal(2,2), vector_normal(2), ident(2,2)

        plane_stress_or_strain = 1
        !! Choose either 1 for plane_stress or 2 for plane_strain

        cebc = 0

        ! loop over all surface loads
        do i = 1, np
          	! if load i is on element e
            if ((loads(i,2) == e) .and. (loads(i,1)==2)) then
              	
				! Load material properties
                thk = mprop(element(e)%mat)%thk
        		young = mprop(element(e)%mat)%young
        		nu = mprop(element(e)%mat)%nu
        		dens = mprop(element(e)%mat)%dens
                
                ! Find coordinates
                nen = element(e)%numnode
                do n = 1, nen
                     xe(2*n-1) = x(element(e)%ix(n),1)
                     xe(2*n  ) = x(element(e)%ix(n),2)
                end do

                ! Find face that the BC is applied to
                eface = loads(i,3)
        		
				! Find surface dimensions
                aa = (xe(3)-xe(1))/2
                bb = (xe(8)-xe(2))/2
        
                nface = 0
                vector_normal = 0
                if (eface == 1) then
                    nface(1,1) = aa
                    nface(1,3) = aa
                    nface(2,2) = aa
                    nface(2,4) = aa
					vector_normal(2) = -1
                elseif (eface == 2) then
                    nface(1,3) = bb
                    nface(1,5) = bb
                    nface(2,4) = bb
                    nface(2,6) = bb
                    vector_normal(1) = 1
                elseif (eface == 3) then
                    nface(1,5) = aa
                    nface(1,7) = aa
                    nface(2,6) = aa
                    nface(2,8) = aa
                    vector_normal(2) = 1
                elseif (eface == 4) then
                    nface(1,1) = bb
                    nface(1,7) = bb
                    nface(2,2) = bb
                    nface(2,8) = bb
                    vector_normal(1) = -1
                endif
				surface_normal(1,1) = vector_normal(1)**2
                surface_normal(1,2) = vector_normal(1)*vector_normal(2)
                surface_normal(2,1) = vector_normal(1)*vector_normal(2)
                surface_normal(2,2) = vector_normal(2)**2
				
				! Calculate damping coefficients
                cs = sqrt(young/(2*(1+nu)*dens))
                select case (plane_stress_or_strain)
                    case (1)
                        ! case of plane stress
                        cp = sqrt(young/((1-nu**2)*dens))
                    case (2)
                        !case of plane strain
                        cp = sqrt(young*(1-nu)/((1+nu)*(1-2*nu)*dens))
                end select
				
				! Build identity matrix
                call identity(ident)

                ! Build z damping matrix
                z = dens*cp*surface_normal + dens*cs*(ident-surface_normal)*loads(i,4) !damping scales with set pressure
                
                ! Build element damping matrix for BCs
                cebc = cebc + matmul(matmul(transpose(nface), z),nface) * thk
                
            end if
        end do
        
    end subroutine plane42_cebc
!
!--------------------------------------------------------------------------------------------------
!
    subroutine plane42_ss(xe, de, young, nu, estress, estrain)

        !! This subrotuine computes the element stress and strain (The location inside the element
        !! where stress and and strain is evaluated, is defined inside the subroutine).

        use fedata

        real(wp), intent(in) :: young
            !! Young's Modulus for this element
        real(wp), intent(in) :: nu 
            !! Poisson's Ratio for this element
        real(wp), dimension(:), intent(in)  :: xe
            !! Nodal coordinates of this element in undeformed configuration (see also [[plane42_ke]])
        real(wp), dimension(:), intent(in)  :: de
            !! Displacement field
            !!
            !! * `de(1:2)` = displacement of element node 1 along \(x\)- and \(y\)-axis, respectively
            !! * `de(3:4)` = displacement of element node 2 along \(x\)- and \(y\)-axis, respectively
            !! * etc...
        real(wp), dimension(:), intent(out) :: estress
            !! Stress at a point inside the element
            !!
            !! * `estress(1)` = \(\sigma_{11}\)
            !! * `estress(2)` = \(\sigma_{22}\)
            !! * `estress(3)` = \(\sigma_{12}\)
        real(wp), dimension(:), intent(out) :: estrain
            !! Strain at a point inside the element
            !!
            !! * `estrain(1)` = \(\epsilon_{11}\)
            !! * `estrain(2)` = \(\epsilon_{22}\)
            !! * `estrain(3)` = \(\epsilon_{12}\)
        real(wp) :: B(3, 8), cmat(3, 3), Ntilde(4,8)
        integer :: i,j
        real(wp) :: loc(3), w(3), xi, eta, det_J, factc

		! Define variable Gauss points and weighting
        loc = 0
        w = 0
        select case (no_gauss_points_s)
            case (1)
            	loc(1) = 0
                w(1) = 2
            case(2)
                loc(1) = 1_wp/sqrt(3.0_wp)
                loc(2) = -1_wp/sqrt(3.0_wp)
                w(1) = 1
                w(2) = 1
            case(3)
            	loc(1) = sqrt(0.6_wp)
                loc(2) = 0
                loc(3) = -sqrt(0.6_wp)
                w(1) = 5.0_wp/9.0_wp
                w(2) = 8.0_wp/9.0_wp
                w(3) = -w(1)
            end select

        ! Calculate strain for each gauss point
        estrain = 0
        do i = 1, no_gauss_points_s
            do j = 1, no_gauss_points_s
                xi = loc(i)
                eta = loc(j)
                ! Calculate B with shape subroutine
                B = 0
				call plane42_shape(xi, eta, Ntilde, B, det_J, xe)
        		! Compute element strain
        		estrain = estrain + matmul(B, de)
        	end do
        end do

        ! Build constitutive matrix (plane stress)
        cmat = 0
        factc = young/(1-nu**2)
        cmat(1,1) = factc
        cmat(1,2) = factc*nu
        cmat(2,1) = factc*nu
        cmat(2,2) = factc
        cmat(3,3) = factc*(1-nu)/2

        ! Compute element stress
        estress = matmul(cmat, estrain)

        ! Compute principal stress and direction
        ! print *, "stress:",estress
    end subroutine plane42_ss

end module plane42
